# Simple MeVD-GRN architecture diagram (SVG), Okabe-Ito colorblind-safe palette. Few words on purpose.
BLUE, ORANGE, GREEN, PURPLE, GREY, INK, BG = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#8A94A0", "#1F2937", "#FFFFFF"
W, H = 1600, 900
o = []
a = o.append
a(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Poppins, Inter, Arial, sans-serif">')
a('<defs><marker id="arr" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#6B7280"/></marker>'
  '<filter id="sh" x="-10%" y="-10%" width="120%" height="130%"><feDropShadow dx="0" dy="2" stdDeviation="3" flood-color="#000" flood-opacity="0.10"/></filter></defs>')
a(f'<rect width="{W}" height="{H}" fill="{BG}"/>')

def card(x, y, w, h, color, title, sub=None, fill=None):
    a(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{fill or "#FFFFFF"}" stroke="{color}" stroke-width="3" filter="url(#sh)"/>')
    a(f'<text x="{x+w/2}" y="{y+h/2+(-4 if sub else 7)}" text-anchor="middle" font-size="22" font-weight="600" fill="{INK}">{title}</text>')
    if sub: a(f'<text x="{x+w/2}" y="{y+h/2+22}" text-anchor="middle" font-size="15" fill="{GREY}">{sub}</text>')
def arrow(x1, y1, x2, y2, dash=False):
    a(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#6B7280" stroke-width="3" marker-end="url(#arr)"' + (' stroke-dasharray="7 6"' if dash else '') + '/>')
def head(x, y, t):
    a(f'<text x="{x}" y="{y}" text-anchor="middle" font-size="17" font-weight="600" letter-spacing="1.5" fill="{GREY}">{t}</text>')

# title
a(f'<text x="60" y="70" font-size="40" font-weight="700" fill="{INK}">MeVD-GRN</text>')
a(f'<text x="60" y="104" font-size="19" fill="{GREY}">Which transcription factor (TF) regulates which gene?</text>')

# column headers
for x, t in [(190, "PAIRED DATA"), (500, "GENE FEATURES"), (820, "GENE GRAPHS"), (1160, "ENCODERS"), (1450, "SCORE")]:
    head(x, 175, t)

# col1: inputs
card(70, 205, 240, 120, BLUE, "scRNA-seq", "expression")
card(70, 385, 240, 120, ORANGE, "scATAC-seq", "open chromatin")
a(f'<path d="M322,215 Q344,365 322,515" fill="none" stroke="{GREY}" stroke-width="3" stroke-linecap="round"/>')
a(f'<text x="196" y="545" text-anchor="middle" font-size="15" fill="{GREY}">same cells</text>')

# col2: features
card(390, 205, 220, 70, BLUE, "Expression")
card(390, 300, 220, 70, ORANGE, "Accessibility")
card(390, 395, 220, 70, PURPLE, "Geneformer", "pretrained prior")
arrow(312, 265, 388, 240); arrow(312, 445, 388, 335)
a(f'<text x="540" y="505" text-anchor="middle" font-size="15" fill="{GREY}">per gene</text>')

# col3: two graphs (mini networks)
def mini(cx, cy, color, label, star=False):
    pts = [(-60,-34),(-12,-52),(46,-30),(-40,22),(16,12),(62,36),(-4,50)]
    edges = [(0,1),(1,2),(0,3),(1,4),(2,4),(3,4),(4,5),(4,6),(3,6)]
    if star: edges = [(4,0),(4,1),(4,2),(4,3),(4,5),(4,6)]
    a(f'<rect x="{cx-110}" y="{cy-95}" width="220" height="180" rx="18" fill="#FFFFFF" stroke="{color}" stroke-width="3" filter="url(#sh)"/>')
    for i, j in edges:
        a(f'<line x1="{cx+pts[i][0]}" y1="{cy-15+pts[i][1]}" x2="{cx+pts[j][0]}" y2="{cy-15+pts[j][1]}" stroke="{color}" stroke-width="2.5" opacity="0.65"/>')
    for k, (px, py) in enumerate(pts):
        r, fc = (11, color) if (star and k == 4) else (7, "#FFFFFF")
        a(f'<circle cx="{cx+px}" cy="{cy-15+py}" r="{r}" fill="{fc}" stroke="{color}" stroke-width="3"/>')
    a(f'<text x="{cx}" y="{cy+72}" text-anchor="middle" font-size="18" font-weight="600" fill="{INK}">{label}</text>')
mini(820, 290, BLUE, "Co-expressed")
mini(820, 520, ORANGE, "TF → likely targets", star=True)
a(f'<text x="820" y="640" text-anchor="middle" font-size="15" fill="{GREY}">built without any labels</text>')
arrow(612, 245, 706, 275); arrow(612, 345, 706, 470)

# col4: towers + embeddings bar
card(1020, 215, 205, 110, BLUE, "RNA tower", "graph network")
card(1020, 395, 205, 110, ORANGE, "ATAC tower", "graph network")
arrow(932, 290, 1018, 270); arrow(932, 505, 1018, 450)
a(f'<rect x="1250" y="215" width="46" height="290" rx="23" fill="{PURPLE}" opacity="0.18" stroke="{PURPLE}" stroke-width="3"/>')
a(f'<text transform="translate(1281,360) rotate(-90)" text-anchor="middle" font-size="18" font-weight="600" fill="{INK}">gene embeddings</text>')
arrow(1227, 270, 1248, 270); arrow(1227, 450, 1248, 450)
# Geneformer prior also feeds the encoders (dashed route underneath)
a(f'<path d="M450,467 L450,672 L1273,672 L1273,508" fill="none" stroke="{PURPLE}" stroke-width="3" stroke-dasharray="7 6" marker-end="url(#arr)"/>')
a(f'<text x="860" y="664" text-anchor="middle" font-size="14" fill="{PURPLE}">also feeds the encoders</text>')

# col5: decoder + output
a(f'<rect x="1340" y="215" width="220" height="394" rx="18" fill="#FFFFFF" stroke="{INK}" stroke-width="3" filter="url(#sh)"/>')
a(f'<circle cx="1393" cy="300" r="26" fill="{INK}"/><text x="1393" y="307" text-anchor="middle" font-size="20" font-weight="700" fill="#FFFFFF">TF</text>')
a(f'<circle cx="1507" cy="300" r="26" fill="#FFFFFF" stroke="{INK}" stroke-width="3"/><text x="1507" y="307" text-anchor="middle" font-size="16" font-weight="600" fill="{INK}">gene</text>')
a(f'<line x1="1421" y1="300" x2="1479" y2="300" stroke="{INK}" stroke-width="3" stroke-dasharray="6 5"/>')
a(f'<rect x="1437" y="268" width="26" height="26" rx="6" fill="{ORANGE}"/>')
a(f'<text x="1450" y="352" text-anchor="middle" font-size="16" fill="{GREY}">gated by openness</text>')
a(f'<text x="1450" y="446" text-anchor="middle" font-size="60" font-weight="700" fill="{GREEN}">0.93</text>')
a(f'<text x="1450" y="478" text-anchor="middle" font-size="17" fill="{INK}">chance TF → gene</text>')
for _y, _t in [(508, "+ hub prior"), (548, "+ motif match")]:
    a(f'<rect x="1368" y="{_y}" width="164" height="32" rx="16" fill="#FFFFFF" stroke="{PURPLE}" stroke-width="2.5" stroke-dasharray="6 4"/>')
    a(f'<text x="1450" y="{_y+22}" text-anchor="middle" font-size="16" font-weight="600" fill="{INK}">{_t}</text>')
a(f'<text x="1450" y="598" text-anchor="middle" font-size="12" fill="{GREY}">add-ons (selected model)</text>')
arrow(1298, 360, 1338, 360)

# bottom: how it learns
a(f'<rect x="60" y="690" width="1500" height="170" rx="22" fill="#F3F4F6"/>')
a(f'<text x="90" y="728" font-size="17" font-weight="600" letter-spacing="1.5" fill="{GREY}">HOW IT LEARNS</text>')
def step(x, color, num, t1, t2):
    a(f'<rect x="{x}" y="745" width="330" height="90" rx="16" fill="#FFFFFF" stroke="{color}" stroke-width="3"/>')
    a(f'<circle cx="{x+40}" cy="790" r="22" fill="{color}"/><text x="{x+40}" y="798" text-anchor="middle" font-size="22" font-weight="700" fill="#FFFFFF">{num}</text>')
    a(f'<text x="{x+78}" y="785" font-size="20" font-weight="600" fill="{INK}">{t1}</text>')
    a(f'<text x="{x+78}" y="810" font-size="15" fill="{GREY}">{t2}</text>')
step(90, BLUE, "1", "ChIP-seq", "TF binds the gene")
step(470, ORANGE, "2", "Knockout", "TF changes the gene")
step(850, GREEN, "3", "ChIP ∩ knockout", "test only: never trained")
arrow(422, 790, 466, 790); arrow(802, 790, 846, 790)
# badges
def badge(x, t, color):
    a(f'<rect x="{x}" y="755" width="290" height="32" rx="16" fill="{color}" opacity="0.16" stroke="{color}" stroke-width="2"/>')
    a(f'<text x="{x+145}" y="777" text-anchor="middle" font-size="16" font-weight="600" fill="{INK}">{t}</text>')
badge(1230, "test pairs never trained on", GREEN)
badge(1230, "", GREEN) if False else None
a(f'<rect x="1230" y="801" width="290" height="32" rx="16" fill="{GREEN}" opacity="0.16" stroke="{GREEN}" stroke-width="2"/>')
a(f'<text x="1375" y="823" text-anchor="middle" font-size="16" font-weight="600" fill="{INK}">zero test → train leakage</text>')
a('</svg>')
open("docs/figures/architecture_simple_2026-10-09.svg", "w").write("\n".join(o))
print("svg written")
