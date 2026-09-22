"""Offline source-provenance regressions for either packaged engine and the R twin.

SCMA_TEST_PACKAGE_ROOT may point to another scmarkeragent package directory.
Run: python -m unittest discover -s tests -p 'test_borrowed_sources.py'
"""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('SCMA_TEST_PACKAGE_ROOT', next(
    (str(p) for p in (BASE / 'src/scmarkeragent', BASE / 'scmarkeragent',
                     BASE / 'pipeline/engine/scmarkeragent') if p.is_dir()), ''
)))
sys.path.insert(0, str(ROOT.parent))
from scmarkeragent import borrowed_context, marker_sources, reporting
from scmarkeragent.ortho_map import OrthoMap


class TestBorrowedSources(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)
        cols = ['species', 'tissue_type', 'disease_normalized', 'cell_type',
                'gene_symbol', 'n_pub_support', 'pmcid', 'pmid', 'source']
        rows = [
            ['Human', 'liver', 'Normal', 'native', 'ALB', 3, 'PMC1', '1', 'native evidence'],
            ['Human', 'pleura', 'Normal', 'borrowed', 'MSLN', 5, 'PMC2', '2', 'cross tissue evidence'],
            ['Mouse', 'heart', 'Normal', 'donated', 'Donor', 7, 'PMC3', '3', 'mouse evidence'],
            ['Rat', 'lung', 'Normal', 'donated', 'Donor', 6, 'PMC4', '4', 'rat evidence'],
            ['Mouse', 'heart', 'Cancer', 'donated', 'Donor', 100, 'PMC5', '5', 'wrong disease'],
            ['Mouse', 'heart', 'Normal', 'donated', 'Ambiguous', 99, 'PMC6', '6', 'ambiguous ortholog'],
            ['Mouse', 'heart', 'Normal', 'donated', 'Unmapped', 99, 'PMC7', '7', 'no ortholog'],
            ['Human', 'pleura', 'Normal', 'unborrowed', 'MSLN', 9, 'PMC8', '8', 'out of scope'],
        ]
        pd.DataFrame(rows, columns=cols).to_csv(self.path / 'sources.csv', index=False)
        for species in ['Mouse', 'Rat']:
            pd.DataFrame([['TARGET', 'Donor'], ['A', 'Ambiguous'], ['B', 'Ambiguous']],
                         columns=['human', 'target']).to_csv(
                             self.path / f'ortho_Human_to_{species}.csv', index=False)
        self.db = marker_sources.SourceDB(str(self.path / 'sources.csv'))
        self.ortho = OrthoMap(str(self.path))
        candidates = []
        for name, gene, borrow in [('native', 'ALB', None),
                                    ('borrowed', 'MSLN', {'donor_species': None}),
                                    ('donated', 'TARGET', {'donor_species': 'Mouse'}),
                                    ('unborrowed', 'MSLN', None)]:
            candidates.append({'cell_type': name, 'claim_role': 'selected',
                'borrowed_context': borrow,
                'decisive_marker_measurements': [{'gene': gene, 'polarity': 'positive',
                    'detection_fraction_in': 0.8, 'detection_fraction_out': 0.1,
                    'avg_log2FC': 2.0, 'auc': 0.9, 'cross_cluster_percentile': 0.75,
                    'publication_support': 7, 'evidence_tier': 'high'}]})
        self.fr = {'res': {'0': {'annotation': 'borrowed', 'candidate_entries': candidates}},
            'scoring': {'context': {'species': 'Human', 'tissue': 'liver', 'disease': 'Normal'},
                        'marker_specificity': {'ALB': 1.0, 'MSLN': 1.0, 'TARGET': 1.0}},
            'dm': {'de': pd.DataFrame({'group': ['0'] * 3,
                'feature': ['ALB', 'MSLN', 'TARGET'], 'padj': [0.001] * 3})}}

    def build(self, fr=None):
        with patch.object(reporting, 'SourceDB', return_value=self.db), \
             patch.object(reporting, 'OrthoMap', return_value=self.ortho), \
             patch.object(marker_sources, 'tissue_members', side_effect=lambda t: [t]), \
             patch.object(marker_sources, 'CROSS_SPECIES', False):
            return reporting.build_marker_evidence('fixture', fr or self.fr)

    def test_export_retains_sources_without_changing_measurements(self):
        after = self.build()
        before_fr = copy.deepcopy(self.fr)
        for entry in before_fr['res']['0']['candidate_entries']:
            entry['borrowed_context'] = None
        before = self.build(before_fr)
        source_cols = ['pmid', 'pmcid', 'source_sentence']
        pd.testing.assert_frame_equal(before.drop(columns=source_cols), after.drop(columns=source_cols))
        indexed = after.set_index('candidate_annotation')
        self.assertEqual(indexed.loc['borrowed', source_cols].tolist(), ['2', 'PMC2', 'cross tissue evidence'])
        self.assertEqual(indexed.loc['donated', source_cols].tolist(), ['3', 'PMC3', 'mouse evidence'])
        self.assertEqual(indexed.loc['native', 'pmcid'], 'PMC1')
        self.assertEqual(indexed.loc['unborrowed', 'pmcid'], 'N/A')
        self.assertEqual(before.set_index('candidate_annotation').loc['borrowed', 'pmcid'], 'N/A')

    def test_donor_disease_and_one_to_one_constraints(self):
        extra = self.db.cell_types_across_tissues('Human', ['donated'], 'Normal',
            donor_species_by_name={'donated': {'Mouse', 'Rat'}}, ortho=self.ortho)
        self.assertEqual({r['pmcid'] for r in extra.all_records('donated', 'TARGET')}, {'PMC3', 'PMC4'})
        for gene in ['A', 'B', 'Unmapped']:
            self.assertEqual(extra.all_records('donated', gene), [])
        native_only = self.db.cell_types_across_tissues('Human', ['donated'], 'Normal')
        self.assertEqual(native_only.all_records('donated', 'TARGET'), [])

    def test_multiple_donors_survive_event_and_result_collection(self):
        events = {str(i): {'borrowed': [{'cell_type': 'donated', 'donor_species': donor}]}
                  for i, donor in enumerate(['Mouse', 'Rat', None, 'Mouse'])}
        self.assertEqual(borrowed_context.borrowed_donor_species(events), {'donated': {'Mouse', 'Rat'}})
        self.fr['res']['1'] = copy.deepcopy(self.fr['res']['0'])
        self.fr['res']['1']['candidate_entries'][2]['borrowed_context']['donor_species'] = 'Rat'
        names, donors = reporting._borrowed_context_names(self.fr)
        self.assertEqual(names, {'borrowed', 'donated'})
        self.assertEqual(donors, {'donated': {'Mouse', 'Rat'}})

    def test_r_report_matches_python(self):
        rscript = shutil.which('Rscript')
        if not rscript:
            self.skipTest('Rscript is not installed')
        payload = copy.deepcopy(self.fr)
        payload['dm']['de'] = payload['dm']['de'].to_dict('records')
        (self.path / 'fr.json').write_text(json.dumps(payload), encoding='utf-8')
        completed = subprocess.run([rscript, str(Path(__file__).with_suffix('.R')),
                                    str(ROOT), str(self.path)], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        r_table = pd.read_csv(self.path / 'r_table.csv', dtype=str, keep_default_na=False)
        pd.testing.assert_frame_equal(self.build().astype(str), r_table)


if __name__ == '__main__':
    unittest.main()
