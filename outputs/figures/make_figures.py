import sys
import os
B=os.path.abspath(os.path.join(os.path.dirname(__file__),"..",".."))
sys.path.insert(0, os.path.join(B,"scripts"))
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from shp_dbf import read_dbf, read_shp_polygons
OUT=os.path.dirname(os.path.abspath(__file__))
_, rows = read_dbf(f"{B}/data/raw/census/county_boundaries/cb_2022_us_county_20m.dbf")
feats = read_shp_polygons(f"{B}/data/raw/census/county_boundaries/cb_2022_us_county_20m.shp")
d = pd.read_csv(f"{B}/outputs/county_year_panel_full.csv", dtype={"county_fips":str})
cty = d.groupby("county_fips").agg(cei=("cei","mean"), ww=("well_water_fraction","first"))
EXCL={"02","15","72","60","66","69","78"}
polys=[];fips=[]
for r,f in zip(rows,feats):
    if r["STATEFP"] in EXCL: continue
    for ring in f:
        polys.append(ring); fips.append(r["GEOID"])
fips=pd.Series(fips)
GREY="#e4e3df"; INK="#2b2b2a"
def cmap(hi): return LinearSegmentedColormap.from_list("s",["#eef3fa",hi])
def draw(ax, vals, cm, vmin, vmax, nodata=GREY):
    v = fips.map(vals).values.astype(float)
    colors=[nodata if np.isnan(x) else cm((x-vmin)/(vmax-vmin)) for x in v]
    pc=PolyCollection(polys, facecolors=colors, edgecolors="white", linewidths=0.15)
    ax.add_collection(pc); ax.set_xlim(-125,-66.5); ax.set_ylim(24,50); ax.set_aspect(1.25); ax.axis("off")
def cbar(fig, ax, cm, vmin, vmax, label):
    sm=plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(vmin,vmax))
    cb=fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.035, pad=0.01, shrink=0.5)
    cb.set_label(label, color=INK, fontsize=9); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=8, colors=INK)
# Fig 1: CEI
fig,ax=plt.subplots(figsize=(8,5.2),dpi=300)
cm=cmap("#1f5fae"); draw(ax,cty.cei.dropna(),cm,0,1); cbar(fig,ax,cm,0,1,"Cumulative Exposure Intensity (mean fraction of observed weeks above WHO 10 µg/L threshold)")
ax.text(-125,24.3,"Grey: no CyAN-matched lake ≥1 km² (excluded from model)",fontsize=8,color="#52514e")
fig.savefig(f"{OUT}/fig1_cei_map.png",bbox_inches="tight",facecolor="white"); plt.close(fig)
# Fig 2: well water
fig,ax=plt.subplots(figsize=(8,5.2),dpi=300)
cm=cmap("#a8470f"); draw(ax,cty.ww.dropna(),cm,0,1); cbar(fig,ax,cm,0,1,"Private/self-supplied well water fraction of population (USGS 2015)")
ax.text(-125,24.3,"Grey: no USGS water-use record",fontsize=8,color="#52514e")
fig.savefig(f"{OUT}/fig2_wellwater_map.png",bbox_inches="tight",facecolor="white"); plt.close(fig)
# Fig 4: forest plot (MI-pooled RR per unit, from pooled_fixed_effects.txt)
terms=[("Bloom exposure (CEI)\nper full 0–1 range",0.997,0.944,1.053),("Well-water fraction\nper 1.0 (0–100%)",1.541,1.392,1.706),
("Poverty rate\nper 1.0 (0–100%)",0.588,0.328,1.055),("Median household income\nper $10,000",0.861,0.846,0.876),("Lake-surface temperature\nper 10 °C",1.008,0.972,1.045)]
fig,ax=plt.subplots(figsize=(8.5,4.3),dpi=400)
n=len(terms)
for i,(nm,rr,lo,hi) in enumerate(terms[::-1]):
    sig = lo>1 or hi<1; c="#1f5fae" if sig else "#6b6a66"
    ax.plot([lo,hi],[i,i],color=c,lw=2.4,solid_capstyle="round"); ax.plot(rr,i,"o",ms=8,color=c,mec="white",mew=1.5)
    ax.annotate(f"{rr:.3f}  ({lo:.3f}–{hi:.3f})",xy=(1.03,i),xycoords=("axes fraction","data"),va="center",fontsize=9,color=INK)
ax.axvline(1,color="#9b9a95",lw=1,ls="--"); ax.set_yticks(range(n)); ax.set_yticklabels([t[0] for t in terms[::-1]],fontsize=9.5,color=INK)
ax.tick_params(axis="y",length=0)
ax.set_xscale("log"); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xticks([0.4,0.6,0.8,1,1.5]); ax.set_xticklabels(["0.4","0.6","0.8","1.0","1.5"],fontsize=9); ax.set_xlim(0.3,1.75); ax.set_ylim(-0.6,n-0.4)
ax.set_xlabel("Rate ratio (95% credible interval), log scale",fontsize=10,color=INK)
ax.annotate("Rate ratio (95% CI)",xy=(1.03,n-0.45),xycoords=("axes fraction","data"),va="bottom",fontsize=9,color="#52514e",fontweight="bold")
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/fig4_forest_plot.png",bbox_inches="tight",facecolor="white"); plt.close(fig)
# Fig 3: suppression
sup=d.groupby("county_fips").agg(s=("deaths_suppressed","mean"),lake=("cei",lambda x:x.notna().any()))
allc=sup.s
from matplotlib.collections import PolyCollection as PC
fig,ax=plt.subplots(figsize=(8,5.2),dpi=300)
cm=cmap("#a8470f")
v=fips.map(allc).values.astype(float)
colors=[GREY if np.isnan(x) else cm(x) for x in v]
lk=fips.map(sup.lake).fillna(False).astype(bool).values
ax.add_collection(PC(polys,facecolors=colors,edgecolors=["#1f5fae" if l else "white" for l in lk],linewidths=[0.45 if l else 0.1 for l in lk]))
ax.set_xlim(-125,-66.5); ax.set_ylim(24,50); ax.set_aspect(1.25); ax.axis("off")
cbar(fig,ax,cm,0,1,"Fraction of 7 study years (2018–2024) with suppressed CDC WONDER mortality count")
ax.text(-125,24.3,"Blue outline: county has CyAN-matched lake/reservoir ≥1 km². Grey: no CDC WONDER rows.",fontsize=8,color="#52514e")
fig.savefig(f"{OUT}/fig3_suppression_map.png",bbox_inches="tight",facecolor="white"); plt.close(fig)
print("done")
