# Simple MeVD-GRN architecture diagram for the scMultiomeGRN comparison (undirected TF-TF link prediction).
# Same visual language as make_architecture_simple.py (Okabe-Ito palette, few words). Writes SVG and PNG.
BLUE, ORANGE, GREEN, PURPLE, GREY, INK, BG = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#8A94A0", "#1F2937", "#FFFFFF"
W, H = 1600, 900
o = []
a = o.append
a(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="DejaVu Sans, Arial, sans-serif">')
a('<defs><marker id="arr" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#6B7280"/></marker>'
  '<marker id="arrs" markerWidth="10" markerHeight="10" refX="2" refY="5" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#1F2937"/></marker></defs>')
a(f'<rect width="{W}" height="{H}" fill="{BG}"/>')


def card(x, y, w, h, color, title, sub=None, dashed=False):
    d = ' stroke-dasharray="8 6"' if dashed else ""
    a(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="#FFFFFF" stroke="{color}" stroke-width="3"{d}/>')
    a(f'<text x="{x+w/2}" y="{y+h/2+(-2 if sub else 8)}" text-anchor="middle" font-size="22" font-weight="700" fill="{INK}">{title}</text>')
    if sub:
        a(f'<text x="{x+w/2}" y="{y+h/2+24}" text-anchor="middle" font-size="15" fill="{GREY}">{sub}</text>')


def arrow(x1, y1, x2, y2, dash=False):
    a(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#6B7280" stroke-width="3" marker-end="url(#arr)"' + (' stroke-dasharray="7 6"' if dash else '') + '/>')


def head(x, y, t):
    a(f'<text x="{x}" y="{y}" text-anchor="middle" font-size="16" font-weight="700" letter-spacing="1.5" fill="{GREY}">{t}</text>')


a(f'<text x="60" y="66" font-size="38" font-weight="700" fill="{INK}">MeVD-GRN</text>')
a(f'<text x="60" y="100" font-size="19" fill="{GREY}">Do two transcription factors (TFs) form a regulatory link?</text>')

for x, t in [(190, "WHAT IT READS"), (560, "TF GRAPHS"), (930, "TWO GRAPH NETWORKS"), (1425, "SCORE A PAIR")]:
    head(x, 160, t)

# column 1: inputs
card(60, 190, 260, 110, BLUE, "scRNA + scATAC", "one summary per TF")
card(60, 330, 260, 110, PURPLE, "Geneformer", "pretrained gene identity")
card(60, 470, 260, 110, GREY, "Matched features", "scMultiomeGRN's own, optional", dashed=True)
# the two extra per-TF sources join one bracket that feeds the encoders (dashed route below)
a(f'<path d="M322,385 L350,385 L350,525 L322,525" fill="none" stroke="{PURPLE}" stroke-width="3"/>')


def mini(cx, cy, color, label, star=False):
    pts = [(-60, -34), (-12, -52), (46, -30), (-40, 22), (16, 12), (62, 36), (-4, 50)]
    edges = [(0, 1), (1, 2), (0, 3), (1, 4), (2, 4), (3, 4), (4, 5), (4, 6), (3, 6)]
    if star:
        edges = [(4, 0), (4, 1), (4, 2), (4, 3), (4, 5), (4, 6), (0, 1)]
    a(f'<rect x="{cx-110}" y="{cy-95}" width="220" height="180" rx="18" fill="#FFFFFF" stroke="{color}" stroke-width="3"/>')
    for i, j in edges:
        a(f'<line x1="{cx+pts[i][0]}" y1="{cy-15+pts[i][1]}" x2="{cx+pts[j][0]}" y2="{cy-15+pts[j][1]}" stroke="{color}" stroke-width="2.5" opacity="0.65"/>')
    for k, (px, py) in enumerate(pts):
        r, fc = (11, color) if (star and k == 4) else (7, "#FFFFFF")
        a(f'<circle cx="{cx+px}" cy="{cy-15+py}" r="{r}" fill="{fc}" stroke="{color}" stroke-width="3"/>')
    a(f'<text x="{cx}" y="{cy+72}" text-anchor="middle" font-size="16" font-weight="700" fill="{INK}">{label}</text>')


# column 2: graphs (built without labels) + optional training links chip
mini(560, 290, BLUE, "Co-expressed TFs")
mini(560, 490, ORANGE, "Likely partner TFs", star=True)
a(f'<text x="560" y="602" text-anchor="middle" font-size="15" fill="{GREY}">built without any labels</text>')
a(f'<rect x="440" y="620" width="240" height="38" rx="19" fill="#FFFFFF" stroke="{GREEN}" stroke-width="2.5" stroke-dasharray="6 4"/>')
a(f'<text x="560" y="645" text-anchor="middle" font-size="16" font-weight="600" fill="{INK}">+ known training links</text>')
arrow(322, 235, 446, 270)
arrow(322, 270, 446, 440)

# column 3: two towers -> one embedding per TF
card(830, 215, 200, 100, BLUE, "RNA network", "graph network")
card(830, 405, 200, 100, ORANGE, "ATAC network", "graph network")
arrow(672, 285, 828, 265); arrow(672, 485, 828, 455)
a(f'<rect x="1075" y="215" width="46" height="290" rx="23" fill="{PURPLE}" opacity="0.18" stroke="{PURPLE}" stroke-width="3"/>')
a(f'<text transform="translate(1106,360) rotate(-90)" text-anchor="middle" font-size="17" font-weight="700" fill="{INK}">one embedding per TF</text>')
arrow(1032, 265, 1073, 285); arrow(1032, 455, 1073, 435)
# the Geneformer identity also feeds the encoders (dashed route underneath)
a(f'<path d="M350,455 L350,690 L1098,690 L1098,508" fill="none" stroke="{PURPLE}" stroke-width="3" stroke-dasharray="7 6" marker-end="url(#arr)"/>')
a(f'<text x="880" y="680" text-anchor="middle" font-size="14" fill="{PURPLE}">also feed the networks</text>')

# column 4: symmetric pair score
a(f'<rect x="1230" y="190" width="300" height="420" rx="18" fill="#FFFFFF" stroke="{INK}" stroke-width="3"/>')
a(f'<circle cx="1305" cy="290" r="30" fill="{INK}"/><text x="1305" y="298" text-anchor="middle" font-size="22" font-weight="700" fill="#FFFFFF">u</text>')
a(f'<circle cx="1455" cy="290" r="30" fill="{INK}"/><text x="1455" y="298" text-anchor="middle" font-size="22" font-weight="700" fill="#FFFFFF">v</text>')
a(f'<line x1="1366" y1="290" x2="1394" y2="290" stroke="{INK}" stroke-width="3" marker-start="url(#arrs)" marker-end="url(#arrs)"/>')
a(f'<text x="1380" y="356" text-anchor="middle" font-size="17" fill="{GREY}">same score either way</text>')
a(f'<text x="1380" y="470" text-anchor="middle" font-size="64" font-weight="700" fill="{GREEN}">0.93</text>')
a(f'<text x="1380" y="506" text-anchor="middle" font-size="17" fill="{INK}">chance of a link</text>')
a(f'<text x="1380" y="532" text-anchor="middle" font-size="14" fill="{GREY}">(example score)</text>')
arrow(1125, 360, 1228, 360)

# bottom: how it is trained and tested
a(f'<rect x="60" y="735" width="1470" height="130" rx="22" fill="#F3F4F6"/>')
a(f'<text x="90" y="772" font-size="16" font-weight="700" letter-spacing="1.5" fill="{GREY}">HOW IT IS TRAINED AND TESTED</text>')


def step(x, color, num, t1, t2, w=410):
    a(f'<rect x="{x}" y="788" width="{w}" height="58" rx="14" fill="#FFFFFF" stroke="{color}" stroke-width="3"/>')
    a(f'<circle cx="{x+32}" cy="817" r="18" fill="{color}"/><text x="{x+32}" y="824" text-anchor="middle" font-size="19" font-weight="700" fill="#FFFFFF">{num}</text>')
    a(f'<text x="{x+62}" y="811" font-size="17" font-weight="700" fill="{INK}">{t1}</text>')
    a(f'<text x="{x+62}" y="832" font-size="14" fill="{GREY}">{t2}</text>')


step(90, BLUE, "1", "Learn from known links", "a part of the TF-TF links, never the test ones")
step(540, ORANGE, "2", "Hold out a test fold", "same official folds as scMultiomeGRN")
step(990, GREEN, "3", "Score with their metric", "AUROC and AUPR on links never trained on")
arrow(500, 817, 538, 817); arrow(950, 817, 988, 817)
a('</svg>')

import cairosvg
svg = "\n".join(o)
open("docs/figures/architecture_undirected_2026-10-17.svg", "w").write(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to="docs/figures/architecture_undirected_2026-10-17.png", output_width=3200)
print("written svg + png")
