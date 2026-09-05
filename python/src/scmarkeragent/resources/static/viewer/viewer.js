
const payload = JSON.parse(document.getElementById("result-data").textContent);
const evidence = new Map(payload.evidence.map(row => [String(row.cluster_id), row]));
const markerMap = new Map();
payload.markers.forEach(row => {
  const key = String(row.cluster_id);
  if (!markerMap.has(key)) markerMap.set(key, []);
  markerMap.get(key).push(row);
});
const make = (tag, cls, text) => {
  const node = document.createElement(tag);
  if (cls) node.className = cls;
  if (text !== undefined && text !== null) node.textContent = String(text);
  return node;
};
const format = (value, digits = 0) => {
  if (value === "" || value === null || value === undefined || Number.isNaN(Number(value))) return "—";
  return Number(value).toLocaleString(undefined, {minimumFractionDigits: digits, maximumFractionDigits: digits});
};
const isNA = (value) => {
  const text = String(value ?? "").trim();
  return !text || text === "N/A";
};
// Machine tokens stay in the downloadable tables; the page says what happened.
// An unmapped value falls through unchanged rather than being hidden.
const PHRASES = {
  resolution_status: {
    resolved: "Single identity",
    mixed: "More than one identity",
    unresolved: "No identity supported"
  },
  annotation_source: {
    cluster_annotation: "Decided by the agent",
    relative_score_fallback: "Top of the evidence ranking",
    no_candidate: "No candidate to choose from"
  },
  llm_status: {
    annotated: "Agent reviewed the evidence",
    enabled: "Agent reviewed the evidence",
    disabled: "Agent step turned off",
    skipped_no_credentials: "Agent step skipped"
  },
  annotation_qc: {
    passed: "Held up on review",
    passed_after_revision: "Revised, then held up",
    demoted_to_parent: "Narrowed name not separable — parent reported",
    failed: "Review did not carry it — best answer shown",
    not_checked: "Not reviewed (no model reachable)"
  }
};
const phrase = (field, value) => {
  if (isNA(value)) return "—";
  const token = String(value).trim();
  return (PHRASES[field] && PHRASES[field][token]) || token;
};
const splitList = (value) => {
  if (isNA(value)) return [];
  return String(value).split(";").map(part => part.trim()).filter(Boolean);
};
const cooccurringGenes = (value) => {
  const grouped = new Map();
  if (isNA(value)) return grouped;
  String(value).split(";").forEach(part => {
    const [gene, identity] = part.split("|").map(piece => (piece || "").trim());
    if (!gene || !identity) return;
    if (!grouped.has(identity)) grouped.set(identity, []);
    grouped.get(identity).push(gene);
  });
  return grouped;
};
const padjText = (value) => {
  if (isNA(value)) return "—";
  const number = Number(value);
  if (!Number.isFinite(number)) return "—";
  if (number !== 0 && Math.abs(number) < 0.001) return number.toExponential(1);
  return number.toFixed(3);
};
const CONF_ORDER = {high: 3, medium: 2, low: 1, "": 0, not_applicable: 0, na: 0};
const confChip = (value) => {
  const token = String(value || "").trim().toLowerCase();
  if (!token || token === "not_applicable" || token === "na") return make("span", "conf conf-na", "n/a");
  const cls = token === "high" ? "conf-high" : token === "medium" ? "conf-medium" : "conf-low";
  return make("span", `conf ${cls}`, token);
};
const resolutionChip = (value) => {
  const token = String(value || "").trim();
  const cls = token.replace(/[^a-z_-]/gi, "");
  return make("span", cls ? `chip res-${cls}` : "chip", phrase("resolution_status", token));
};
const clusterBody = document.getElementById("cluster-body");
const query = document.getElementById("cluster-search");
const resolution = document.getElementById("resolution-filter");
const count = document.getElementById("cluster-count");
const detail = document.getElementById("cluster-detail");
let sortState = {key: "cluster_id", ascending: true};
let selectedCluster = null;
let detailFilters = {candidate: "", polarity: "positive", role: "all"};

function addSummaryCards() {
  const host = document.querySelector(".cards");
  if (!host) return;
  const clusters = payload.clusters;
  const noCandidate = clusters.filter(row => String(row.annotation_source || "").trim() === "no_candidate").length;
  const identities = new Set(
    clusters.map(row => String(row.primary_annotation || "").trim()).filter(name => name && name !== "Unknown" && name !== "N/A")
  ).size;
  const mixed = clusters.filter(row => String(row.resolution_status || "").trim() === "mixed").length;
  const abstained = clusters.filter(row => String(row.primary_annotation || "").trim() === "Unknown").length;
  const flagged = clusters.filter(row => !isNA(row.claim_warnings)).length;
  [
    ["With candidates", `${clusters.length - noCandidate} / ${clusters.length}`, "clusters with an identity to choose from"],
    ["Cell types", identities.toLocaleString(), "distinct identities called"],
    ["Mixed", mixed.toLocaleString(), "carrying more than one identity"],
    ["Abstained", abstained.toLocaleString(), "left as Unknown"],
    ["Audit flags", flagged.toLocaleString(), "clusters with claim warnings"]
  ].forEach(([label, value, note]) => {
    const card = make("div", "card");
    card.append(make("div", "metric-label", label), make("div", "metric-value", value), make("div", "metric-note", note));
    host.append(card);
  });
}

function addDefinition(list, label, value) {
  list.append(make("dt", "", label), make("dd", "", value || "—"));
}
function candidateEntries(records) {
  const seen = new Map();
  records.forEach(row => {
    const name = String(row.candidate_annotation || "").trim();
    if (!name) return;
    if (!seen.has(name)) {
      seen.set(name, {
        name,
        selected: String(row.is_selected_annotation).toLowerCase() === "true",
        count: 0
      });
    }
    seen.get(name).count += 1;
  });
  return [...seen.values()].sort(
    (a, b) => Number(b.selected) - Number(a.selected) || a.name.localeCompare(b.name)
  );
}
function filteredMarkers(records) {
  let list = records;
  if (detailFilters.candidate) {
    list = list.filter(row => String(row.candidate_annotation || "").trim() === detailFilters.candidate);
  }
  if (detailFilters.polarity !== "all") {
    list = list.filter(row => String(row.marker_polarity || "positive") === detailFilters.polarity);
  }
  if (detailFilters.role !== "all") {
    list = list.filter(row => String(row.marker_role || "") === detailFilters.role);
  }
  return list;
}
function roleCountsFor(records) {
  const counts = {decisive: 0, support: 0, panel: 0};
  records
    .filter(row => !detailFilters.candidate || String(row.candidate_annotation || "").trim() === detailFilters.candidate)
    .forEach(row => {
      if (counts[row.marker_role] !== undefined) counts[row.marker_role] += 1;
    });
  return counts;
}
function markerTable(records, sourcesOnly = false) {
  if (!records.length) return make("p", "empty", "No marker rows match this filter.");
  const wrap = make("div", "table-wrap");
  const table = document.createElement("table");
  const head = document.createElement("thead");
  const header = document.createElement("tr");
  const labels = sourcesOnly
    ? ["Gene", "Evidence", "PMCID", "Source sentence"]
    : ["Gene", "Role", "In cluster", "Outside cluster", "avg log2FC", "adj. p-value", "Percentile", "Specificity", "Publications", "Tier", "Evidence"];
  labels.forEach(label => header.append(make("th", "", label)));
  head.append(header); table.append(head);
  const body = document.createElement("tbody");
  records.forEach(row => {
    const tr = document.createElement("tr");
    const geneCell = make("td", "");
    geneCell.append(make("span", "", row.gene));
    // A curated exclusion the cluster detects reads as support unless it is marked as
    // what it is: evidence the assigned identity had to answer, not evidence for it.
    if (row.marker_polarity === "negative") {
      geneCell.append(make("span", "badge exclusion", "exclusion detected"));
    }
    tr.append(geneCell);
    if (sourcesOnly) {
      const evidenceCell = document.createElement("td");
      evidenceCell.append(make("span", "badge", "LLM-validated · DB-cited"));
      tr.append(evidenceCell);
      const pmcidCell = document.createElement("td");
      const link = make("a", "", row.pmcid || "—");
      if (row.pmcid) {
        link.href = `https://www.ncbi.nlm.nih.gov/pmc/articles/${encodeURIComponent(row.pmcid)}/`;
        link.target = "_blank"; link.rel = "noopener";
      }
      pmcidCell.append(link); tr.append(pmcidCell);
      const sourceCell = make("td", "source-cell");
      const details = document.createElement("details");
      details.append(make("summary", "", "Read cited source"));
      details.append(make("p", "", row.source_sentence || "No source sentence available."));
      sourceCell.append(details); tr.append(sourceCell);
    } else {
      const roleCell = document.createElement("td");
      const roleToken = String(row.marker_role || "").trim().replace(/[^a-z_-]/gi, "");
      roleCell.append(make("span", roleToken ? `role role-${roleToken}` : "role", row.marker_role));
      tr.append(roleCell);
      tr.append(make("td", "num", format(row.cluster_detection_fraction, 4)));
      tr.append(make("td", "num", format(row.out_of_cluster_detection_fraction, 4)));
      tr.append(make("td", "num", format(row.average_log2_fold_change, 4)));
      tr.append(make("td", "num", padjText(row.adjusted_p_value)));
      tr.append(make("td", "num", format(row.cross_cluster_percentile, 2)));
      tr.append(make("td", "num", format(row.marker_specificity_weight, 4)));
      tr.append(make("td", "num", format(row.publication_support_count)));
      tr.append(make("td", "", isNA(row.evidence_tier) ? "—" : row.evidence_tier));
      const badgeCell = document.createElement("td");
      badgeCell.append(make("span", "badge", "LLM-validated · DB-cited"));
      tr.append(badgeCell);
    }
    body.append(tr);
  });
  table.append(body); wrap.append(table); return wrap;
}
function evidenceToolbar(records, rerender) {
  const toolbar = make("div", "detail-toolbar");
  const entries = candidateEntries(records);
  if (entries.length > 1) {
    const candidates = make("div", "candidates");
    candidates.append(make("span", "seg-label", "Identity"));
    entries.forEach(entry => {
      const button = make("button", `candidate${entry.name === detailFilters.candidate ? " active" : ""}`);
      button.type = "button";
      button.append(make("span", "", entry.name));
      if (entry.selected) button.append(make("span", "candidate-tag", "assigned"));
      button.append(make("span", "candidate-count", entry.count));
      button.addEventListener("click", () => {
        detailFilters.candidate = entry.name === detailFilters.candidate ? "" : entry.name;
        rerender();
      });
      candidates.append(button);
    });
    toolbar.append(candidates);
  }
  const polarity = make("div", "segmented");
  ["positive", "negative", "all"].forEach(option => {
    const button = make("button", option === detailFilters.polarity ? "active" : "", option);
    button.type = "button";
    button.addEventListener("click", () => {
      detailFilters.polarity = option;
      rerender();
    });
    polarity.append(button);
  });
  const counts = roleCountsFor(records);
  const roles = make("div", "segmented");
  const allButton = make("button", detailFilters.role === "all" ? "active" : "", "all roles");
  allButton.type = "button";
  allButton.addEventListener("click", () => {
    detailFilters.role = "all";
    rerender();
  });
  roles.append(allButton);
  ["decisive", "support", "panel"].forEach(role => {
    const button = make("button", role === detailFilters.role ? "active" : "", `${role} (${counts[role]})`);
    button.type = "button";
    button.disabled = !counts[role];
    button.addEventListener("click", () => {
      detailFilters.role = role;
      rerender();
    });
    roles.append(button);
  });
  toolbar.append(polarity, roles);
  return toolbar;
}
function showDetail(row, scroll = false) {
  selectedCluster = String(row.cluster_id);
  detailFilters = {candidate: "", polarity: "positive", role: "all"};
  const record = evidence.get(selectedCluster);
  const markers = markerMap.get(selectedCluster) || [];
  const assigned = candidateEntries(markers).find(entry => entry.selected);
  if (assigned) detailFilters.candidate = assigned.name;
  document.querySelectorAll("#cluster-body tr").forEach(tr => tr.classList.toggle("selected", tr.dataset.cluster === selectedCluster));
  detail.replaceChildren();
  const header = make("div", "detail-header");
  const heading = make("div");
  heading.append(make("h2", "", `Cluster ${row.cluster_id}`));
  heading.append(make("span", "chip", record?.annotation?.cell_type_annotation || row.cell_type_annotation));
  header.append(heading);
  header.append(make("span", "count", `${markers.length} source-backed marker rows`));
  detail.append(header);
  const tabs = make("div", "tabs detail-tabs");
  [["overview", "Overview"], ["markers", "Identity markers"], ["sources", "Sources"]].forEach(([id, label], index) => {
    const button = make("button", `tab detail-tab${index === 0 ? " active" : ""}`, label);
    button.dataset.view = id; tabs.append(button);
  });
  detail.append(tabs);
  const overview = make("section", "detail-view active"); overview.id = "detail-overview";
  const grid = make("div", "overview-grid");
  const definitions = make("dl", "kv");
  addDefinition(definitions, "Cells", format(row.n_cells));
  addDefinition(definitions, "Cell type", row.cell_type_annotation);
  addDefinition(definitions, "Cell subtype", row.cell_subtype_annotation);
  addDefinition(definitions, "Cell Ontology", row.cell_ontology);
  definitions.append(make("dt", "", "Confidence"));
  const confDd = make("dd"); confDd.append(confChip(row.annotation_confidence)); definitions.append(confDd);
  definitions.append(make("dt", "", "Resolution"));
  const resDd = make("dd"); resDd.append(resolutionChip(row.resolution_status)); definitions.append(resDd);
  addDefinition(definitions, "Decided by", phrase("annotation_source", row.annotation_source));
  addDefinition(definitions, "Agent", phrase("llm_status", row.llm_status));
  addDefinition(definitions, "Quality control", phrase("annotation_qc", row.annotation_qc));
  addDefinition(definitions, "Candidate pool", format(record?.audit?.candidate_pool_size));
  grid.append(definitions);
  const rationale = make("div");
  rationale.append(make("h3", "", "Annotation rationale"));
  rationale.append(make("p", "", record?.evidence?.annotation_rationale || "No rationale available."));
  const groups = record?.evidence?.identity_groups || [];
  if (groups.length) {
    rationale.append(make("h3", "", "Candidate names grouped by identity"));
    groups.forEach(group => rationale.append(make("p", "", `${group.identity}: ${(group.candidates || []).join(", ")}`)));
  }
  if (row.alternative_candidates && row.alternative_candidates !== "N/A") {
    rationale.append(make("h3", "", "Other identities this cluster also carries"));
    rationale.append(make("p", "", row.alternative_candidates));
  }
  if (row.cooccurring_markers && row.cooccurring_markers !== "N/A") {
    rationale.append(make("h3", "", "Markers behind those co-occurring identities"));
    rationale.append(make("p", "", row.cooccurring_markers));
  }
  if (row.pmcid && row.pmcid !== "N/A") {
    rationale.append(make("h3", "", "Publications behind the key markers"));
    rationale.append(make("p", "", row.pmcid));
  }
  if (row.cooccurring_pmcid && row.cooccurring_pmcid !== "N/A") {
    rationale.append(make("h3", "", "Publications behind the co-occurring markers"));
    rationale.append(make("p", "", row.cooccurring_pmcid));
  }
  grid.append(rationale); overview.append(grid); detail.append(overview);
  const markerView = make("section", "detail-view"); markerView.id = "detail-markers";
  const sourceView = make("section", "detail-view"); sourceView.id = "detail-sources";
  if (markers.length) {
    const renderEvidence = () => {
      const rows = filteredMarkers(markers);
      markerView.replaceChildren(evidenceToolbar(markers, renderEvidence), markerTable(rows));
      sourceView.replaceChildren(evidenceToolbar(markers, renderEvidence), markerTable(rows, true));
    };
    renderEvidence();
  } else {
    markerView.append(make("p", "empty", "No LLM-validated, source-backed identity markers were available for this cluster."));
    sourceView.append(make("p", "empty", "No source-backed identity markers were available."));
  }
  detail.append(markerView);
  detail.append(sourceView);
  detail.querySelectorAll(".detail-tab").forEach(button => {
    button.addEventListener("click", () => {
      detail.querySelectorAll(".detail-tab").forEach(item => item.classList.remove("active"));
      detail.querySelectorAll(".detail-view").forEach(item => item.classList.remove("active"));
      button.classList.add("active");
      document.getElementById(`detail-${button.dataset.view}`).classList.add("active");
    });
  });
  if (scroll) detail.scrollIntoView({behavior: "smooth", block: "start"});
}
function renderClusters() {
  const needle = query.value.trim().toLowerCase();
  const selectedResolution = resolution.value;
  const filtered = payload.clusters.filter(row => {
    const text = [row.cluster_id, row.cell_type_annotation, row.cell_subtype_annotation, row.cell_state, row.cell_ontology, row.key_markers, row.cooccurring_markers, row.alternative_candidates].join(" ").toLowerCase();
    return (!needle || text.includes(needle)) && (!selectedResolution || row.resolution_status === selectedResolution);
  });
  filtered.sort((a, b) => {
    let left = a[sortState.key], right = b[sortState.key];
    if (["cluster_id", "n_cells"].includes(sortState.key)) { left = Number(left); right = Number(right); }
    else if (["annotation_confidence"].includes(sortState.key)) {
      left = CONF_ORDER[String(left || "").toLowerCase()] ?? 0;
      right = CONF_ORDER[String(right || "").toLowerCase()] ?? 0;
    }
    const result = left < right ? -1 : left > right ? 1 : 0;
    return sortState.ascending ? result : -result;
  });
  clusterBody.replaceChildren();
  filtered.forEach(row => {
    const tr = document.createElement("tr"); tr.dataset.clickable = "true"; tr.dataset.cluster = String(row.cluster_id);
    tr.append(make("td", "num", row.cluster_id), make("td", "num", format(row.n_cells)));
    const typeCell = document.createElement("td");
    typeCell.append(make("div", "cell-primary", row.cell_type_annotation));
    if (!isNA(row.cell_state)) typeCell.append(make("div", "sub", `state: ${row.cell_state}`));
    const coGenes = cooccurringGenes(row.cooccurring_markers);
    splitList(row.alternative_candidates).forEach(name => {
      const genes = coGenes.get(name) || [];
      let text = `also: ${name}`;
      if (genes.length) {
        text += ` — ${genes.slice(0, 4).join(", ")}`;
        if (genes.length > 4) text += ` +${genes.length - 4} more`;
      }
      typeCell.append(make("div", "sub sub-alt", text));
    });
    if (!isNA(row.unlisted_identity)) {
      typeCell.append(make("div", "sub", `shortlist may have missed: ${row.unlisted_identity}`));
    }
    tr.append(typeCell);
    tr.append(make("td", "", row.cell_subtype_annotation), make("td", "", row.cell_ontology));
    const confCell = document.createElement("td"); confCell.append(confChip(row.annotation_confidence)); tr.append(confCell);
    const state = document.createElement("td");
    state.append(resolutionChip(row.resolution_status));
    if (!isNA(row.claim_warnings)) state.append(make("span", "flag-audit", "audit"));
    const qcToken = String(row.annotation_qc || "").trim();
    if (qcToken && qcToken !== "passed" && !isNA(qcToken)) {
      state.append(make("div", "sub", phrase("annotation_qc", qcToken)));
    }
    tr.append(state);
    tr.append(make("td", "", phrase("annotation_source", row.annotation_source)));
    tr.addEventListener("click", () => showDetail(row, true)); clusterBody.append(tr);
  });
  count.textContent = `${filtered.length} of ${payload.clusters.length} clusters`;
  if (selectedCluster) document.querySelectorAll("#cluster-body tr").forEach(tr => tr.classList.toggle("selected", tr.dataset.cluster === selectedCluster));
}
document.querySelectorAll("[data-sort]").forEach(button => button.addEventListener("click", () => {
  const key = button.dataset.sort;
  sortState = {key, ascending: sortState.key === key ? !sortState.ascending : true};
  renderClusters();
}));
query.addEventListener("input", renderClusters); resolution.addEventListener("change", renderClusters);
document.querySelectorAll("#resolution-filter option").forEach(option => {
  if (option.value) option.textContent = phrase("resolution_status", option.value);
});
document.querySelectorAll(".figure-tab").forEach(button => button.addEventListener("click", () => {
  document.querySelectorAll(".figure-tab").forEach(item => item.classList.remove("active"));
  document.querySelectorAll(".figure-stage img").forEach(item => item.classList.remove("active"));
  button.classList.add("active"); document.getElementById(button.dataset.figure).classList.add("active");
}));
addSummaryCards();
renderClusters();
if (payload.clusters.length) showDetail(payload.clusters[0]);
