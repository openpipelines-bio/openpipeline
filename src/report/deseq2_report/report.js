(function(){
  "use strict";
  var D = window.REPORT_DATA;
  var $ = function(id){ return document.getElementById(id); };
  var fmt = function(n){ return n.toLocaleString("en-GB"); };
  var signed = function(x, d){ return (x > 0 ? "+" : "") + x.toFixed(d); };

  /* ---- theme ---- */
  function currentTheme(){
    var attr = document.documentElement.getAttribute("data-theme");
    if(attr) return attr;
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  function css(v){ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }

  var k = D.kpis;
  var PADJ = (k && k.padj_threshold) || 0.05;
  var LFC = (k && k.log2fc_threshold) || 0;
  var sigRule = "padj &lt; " + PADJ + (LFC ? " and |log2FC| &gt; " + LFC : "");

  /* ---- header / meta / KPIs ---- */
  var joinSet = function(parts){ return parts.filter(function(p){ return p; }).join("  \u00b7  "); };
  if(D.meta.title){ $("title").textContent = D.meta.title; document.title = D.meta.title; }
  $("subtitle").textContent = joinSet([D.meta.project, D.meta.design]);
  $("meta").innerHTML = [["Pipeline", D.meta.pipeline], ["Reference", D.meta.ref], ["Design", D.meta.design]]
    .filter(function(m){ return m[1]; })
    .map(function(m){ return "<span><b>" + m[0] + "</b> " + m[1] + "</span>"; }).join("");
  $("foot-meta").textContent = joinSet([D.meta.pipeline, D.meta.ref, D.meta.attribution]);

  $("lead").innerHTML = "Across <b>" + k.samples + " samples</b>, <b>" + fmt(k.genes_quantified) +
    "</b> genes were quantified. At " + sigRule + ", <b>" + fmt(k.de_genes) +
    "</b> genes are differentially expressed (" + fmt(k.up) + " up, " + fmt(k.down) +
    " down)." + (k.median_mapping!=null ? " Median mapping rate was " + k.median_mapping + "%." : "");

  var tiles = D.tiles || [
    {l:"Samples", v:fmt(k.samples)},
    {l:"Genes quantified", v:fmt(k.genes_quantified)},
    {l:"Median mapping", v:k.median_mapping+"%"},
    {l:"Reads assigned", v:k.assigned+"%"},
    {l:"Duplication", v:k.median_dup+"%"},
    {l:"DE genes (padj<"+PADJ+")", v:fmt(k.de_genes),
     sub:'<span class="pos">&#9650; '+k.up+'</span> &nbsp; <span class="neg">&#9660; '+k.down+'</span>', accent:true}
  ];
  $("kpis").innerHTML = tiles.map(function(t){
    return '<div class="kpi'+(t.accent?' accent':'')+'">' +
      '<div class="k-label">'+t.l+'</div><div class="k-val">'+t.v+'</div>' +
      (t.sub?'<div class="k-sub">'+t.sub+'</div>':'') + '</div>';
  }).join("");

  /* ---- Plotly shared layout ---- */
  var CFG = {displayModeBar:false, responsive:true};
  /* categorical slots for the condition groups; the first three pass all-pairs CVD checks */
  var PAL=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#4a3aa7","#e34948"];
  function axis(extra){
    var a={gridcolor:css("--chart-grid"), zeroline:false, linecolor:css("--chart-grid"), ticks:"outside", tickcolor:css("--chart-grid")};
    for(var key in (extra||{})) a[key]=extra[key];
    return a;
  }
  function baseLayout(extra){
    var paper=css("--chart-paper");
    var lay = {
      paper_bgcolor:paper, plot_bgcolor:paper,
      font:{color:css("--chart-ink"), family:'system-ui,-apple-system,"Segoe UI",sans-serif', size:12},
      margin:{l:52,r:16,t:10,b:44}, hovermode:"closest",
      xaxis:axis(), yaxis:axis(), showlegend:false
    };
    for(var key in (extra||{})) lay[key]=extra[key];
    return lay;
  }
  function groupNames(list){ var names=[]; list.forEach(function(g){ if(names.indexOf(g)<0) names.push(g); }); return names; }
  function groupColor(names, g){ return PAL[names.indexOf(g) % PAL.length]; }

  /* ---- PCA ---- */
  function drawPCA(){
    var p=D.pca, names=groupNames(p.group);
    if(p.hint) $("pca-hint").innerHTML = p.hint;
    var traces=names.map(function(g){
      var xs=[],ys=[],tx=[];
      for(var i=0;i<p.group.length;i++) if(p.group[i]===g){xs.push(p.pc1[i]);ys.push(p.pc2[i]);tx.push(p.sample[i]);}
      return {x:xs,y:ys,text:tx,name:g,mode:"markers",type:"scatter",
        marker:{size:12,color:groupColor(names,g),line:{width:1.5,color:css("--chart-paper")}},
        hovertemplate:"<b>%{text}</b><br>"+g+"<br>PC1 %{x}, PC2 %{y}<extra></extra>"};
    });
    var lay=baseLayout({showlegend:true, legend:{orientation:"h",y:1.12,x:0,font:{size:11}},
      xaxis:axis({title:{text:"PC1 ("+p.var[0]+"%)",font:{size:11}}}),
      yaxis:axis({title:{text:"PC2 ("+p.var[1]+"%)",font:{size:11}}})});
    Plotly.newPlot("pca",traces,lay,CFG);
  }

  /* ---- Volcano and MA share the same points and classes ---- */
  function deTraces(xs, ys, yLabel){
    var v=D.volcano, map={ns:css("--ns"),up:css("--up"),down:css("--down")};
    var groups={ns:[],up:[],down:[]};
    for(var i=0;i<v.x.length;i++) groups[v.cls[i]].push(i);
    function tr(key,name){
      var idx=groups[key];
      return {x:idx.map(function(i){return xs[i];}), y:idx.map(function(i){return ys[i];}),
        text:idx.map(function(i){return v.gene[i];}), name:name, mode:"markers", type:"scatter",
        marker:{size:key==="ns"?4:6, color:map[key], opacity:key==="ns"?0.5:0.85,
          line:{width:key==="ns"?0:0.5, color:css("--chart-paper")}},
        hovertemplate:yLabel};
    }
    return [tr("ns","n.s."),tr("down","Down"),tr("up","Up")];
  }
  function drawVolcano(){
    var v=D.volcano;
    $("volcano-hint").innerHTML = "Hover a point for the gene. Red up, blue down at " + sigRule +
      (k.shrunk ? " (shrunken fold changes)" : "") + ".";
    var lay=baseLayout({
      xaxis:axis({title:{text:"log2 fold change",font:{size:11}},zeroline:true,zerolinecolor:css("--chart-grid")}),
      yaxis:axis({title:{text:"-log10 p-value",font:{size:11}}})});
    Plotly.newPlot("volcano",deTraces(v.x, v.y, "<b>%{text}</b><br>log2FC %{x}<br>-log10 p %{y}<extra></extra>"),lay,CFG);
  }
  function drawMA(){
    var v=D.volcano;
    if(!v.base_mean) return;
    $("ma-card").hidden=false;
    var lay=baseLayout({
      xaxis:axis({title:{text:"mean of normalised counts",font:{size:11}},type:"log"}),
      yaxis:axis({title:{text:"log2 fold change",font:{size:11}},zeroline:true,zerolinecolor:css("--chart-ink")})});
    Plotly.newPlot("ma",deTraces(v.base_mean, v.x, "<b>%{text}</b><br>mean %{x:.0f}<br>log2FC %{y}<extra></extra>"),lay,CFG);
  }

  /* ---- Correlation heatmap ---- */
  function drawCorr(){
    var c=D.corr;
    var scale=[[0,css("--chart-paper")],[0.4,"#bcd6f4"],[0.7,"#5a97e0"],[1,"#1c5cab"]];
    /* scale from the data: floor at the lowest off-diagonal value so an outlier sample stays visible */
    var lo=1; c.z.forEach(function(r,i){ r.forEach(function(v,j){ if(i!==j && v<lo) lo=v; }); });
    var zmin=Math.max(0, Math.floor(lo*20)/20 - 0.05);
    var tr={z:c.z,x:c.labels,y:c.labels,type:"heatmap",colorscale:scale,zmin:zmin,zmax:1,
      xgap:1,ygap:1,colorbar:{thickness:10,len:0.85,tickfont:{size:10},outlinewidth:0},
      hovertemplate:"%{y} \u00d7 %{x}<br>r = %{z}<extra></extra>"};
    var lay=baseLayout({margin:{l:70,r:16,t:10,b:70},
      xaxis:{tickfont:{size:9},tickangle:-45,showgrid:false,linecolor:css("--chart-grid")},
      yaxis:{tickfont:{size:9},autorange:"reversed",showgrid:false,linecolor:css("--chart-grid")}});
    Plotly.newPlot("corr",[tr],lay,CFG);
  }

  /* ---- Top-gene heatmap (diverging z-score); first gene on top, every gene labelled ---- */
  function drawHeat(){
    var h=D.heatmap;
    if(!h.genes.length){
      $("heat-hint").innerHTML = "No significant genes at " + sigRule + ".";
      $("heat").hidden = true;
      return;
    }
    if(D.summary) $("heat-hint").textContent = "Row-scaled (z-score) expression of the genes with the lowest padj, up then down.";
    var rows=h.genes.map(function(_,i){ return i; });
    var text=h.genes.map(function(g){ return h.samples.map(function(){ return g; }); });
    var scale=[[0,"#2a78d6"],[0.5,css("--chart-grid")],[1,"#e34948"]];
    var tr={z:h.z,x:h.samples,y:rows,text:text,type:"heatmap",colorscale:scale,zmin:-2,zmax:2,
      xgap:1,ygap:1,colorbar:{thickness:10,len:0.85,tickfont:{size:10},outlinewidth:0,title:{text:"z",font:{size:10}}},
      hovertemplate:"%{text} in %{x}<br>z = %{z}<extra></extra>"};
    var lay=baseLayout({margin:{l:74,r:16,t:10,b:70},
      xaxis:{tickfont:{size:9},tickangle:-45,showgrid:false,linecolor:css("--chart-grid")},
      yaxis:{tickfont:{size:9},showgrid:false,linecolor:css("--chart-grid"),autorange:"reversed",
        tickmode:"array",tickvals:rows,ticktext:h.genes}});
    Plotly.newPlot("heat",[tr],lay,CFG);
  }

  /* ---- Key genes: small multiples, one per gene, paired samples joined ---- */
  function drawGenes(){
    var G=D.genes;
    if(!G || !G.items || !G.items.length) return;
    $("genes-card").hidden=false;
    $("genes-hint").innerHTML = "Normalised counts per sample (log scale)" +
      (G.pairs ? "; lines join samples from the same " + G.pair_label + "." : ".");
    var names=groupNames(G.groups), n=G.items.length, cols=n>1?2:1, rowsN=Math.ceil(n/cols);
    var traces=[], ann=[], lay=baseLayout({
      grid:{rows:rowsN, columns:cols, pattern:"independent", xgap:0.16, ygap:0.32},
      margin:{l:56,r:12,t:28,b:30}, showlegend:true, legend:{orientation:"h",y:-0.08,x:0,font:{size:11}}});
    G.items.forEach(function(it, gi){
      var sfx = gi===0 ? "" : String(gi+1);
      var xa="x"+sfx, ya="y"+sfx;
      if(G.pairs){
        groupNames(G.pairs).forEach(function(p){
          var xs=[], ys=[];
          for(var i=0;i<G.pairs.length;i++) if(G.pairs[i]===p){ xs.push(G.groups[i]); ys.push(it.values[i]); }
          traces.push({x:xs,y:ys,xaxis:xa,yaxis:ya,mode:"lines",type:"scatter",hoverinfo:"skip",
            line:{width:1.5,color:css("--ns")},showlegend:false});
        });
      }
      names.forEach(function(g){
        var xs=[], ys=[], tx=[];
        for(var i=0;i<G.groups.length;i++) if(G.groups[i]===g){ xs.push(g); ys.push(it.values[i]); tx.push(G.samples[i]); }
        traces.push({x:xs,y:ys,text:tx,xaxis:xa,yaxis:ya,name:g,legendgroup:g,showlegend:gi===0,
          mode:"markers",type:"scatter",
          marker:{size:10,color:groupColor(names,g),line:{width:1.5,color:css("--chart-paper")}},
          hovertemplate:"<b>"+it.gene+"</b> in %{text}<br>%{y:.0f} normalised counts<extra></extra>"});
      });
      lay["xaxis"+sfx]=axis({categoryorder:"array",categoryarray:names,range:[-0.4,names.length-0.6],tickfont:{size:10},ticks:""});
      lay["yaxis"+sfx]=axis({type:"log",dtick:1,tickformat:"~s",tickfont:{size:10},
        title:(gi%cols===0)?{text:"normalised counts",font:{size:10}}:undefined});
      ann.push({xref:xa+" domain",yref:ya+" domain",x:0,y:1.04,xanchor:"left",yanchor:"bottom",showarrow:false,
        text:"<b>"+it.gene+"</b>  log2FC "+signed(it.log2fc,2)+"  \u00b7  padj "+it.padj.toExponential(1),
        font:{size:11,color:css("--ink")}});
    });
    lay.annotations=ann;
    Plotly.newPlot("genes",traces,lay,CFG);
  }

  /* ---- P-value histogram with the estimated null level ---- */
  function drawPvals(){
    var P=D.pvalues;
    if(!P) return;
    $("pval-card").hidden=false;
    var w=P.bin_width, xs=P.counts.map(function(_,i){ return (i+0.5)*w; });
    var nullLevel=P.pi0*P.n*w;
    var nonnull=P.n_nonnull!=null ? P.n_nonnull : Math.round((1-P.pi0)*P.n);
    $("pval-hint").innerHTML = "Raw p-values of the " + fmt(P.n) + " tested genes. A flat floor means the test is well calibrated; " +
      "the spike near 0 is the real signal: about <b>" + fmt(nonnull) + "</b> genes (" + Math.round((1-P.pi0)*100) +
      "%) are estimated to respond.";
    var tr={x:xs,y:P.counts,type:"bar",width:w*0.92,marker:{color:css("--accent")},
      hovertemplate:"p %{customdata}<br>%{y} genes<extra></extra>",
      customdata:xs.map(function(x){ return (x-w/2).toFixed(3)+"\u2013"+(x+w/2).toFixed(3); })};
    var lay=baseLayout({margin:{l:52,r:16,t:10,b:40},
      xaxis:axis({title:{text:"p-value",font:{size:11}},range:[0,1]}),
      yaxis:axis({title:{text:"genes",font:{size:11}}}),
      shapes:[{type:"line",xref:"x",yref:"y",x0:0,x1:1,y0:nullLevel,y1:nullLevel,
        line:{color:css("--chart-ink"),width:1.5,dash:"dash"}}],
      annotations:[{x:1,y:nullLevel,xref:"x",yref:"y",xanchor:"right",yanchor:"bottom",showarrow:false,
        text:"expected if no gene responded (\u03c0\u2080 = "+P.pi0.toFixed(2)+")",font:{size:11,color:css("--chart-ink")}}]});
    Plotly.newPlot("pval",[tr],lay,CFG);
  }

  /* ---- results table ---- */
  var rows = D.table.slice();
  var sortK="padj", sortDir=1, filter="";
  function dirOf(r){ return (r.padj >= PADJ || Math.abs(r.log2fc) <= LFC) ? "ns" : (r.log2fc>0 ? "up":"down"); }
  function render(){
    var f=filter.trim().toUpperCase();
    var view=rows.filter(function(r){ return !f || r.gene.toUpperCase().indexOf(f)>=0; });
    view.sort(function(a,b){
      var av,bv;
      if(sortK==="gene"){av=a.gene;bv=b.gene;return av<bv?-sortDir:av>bv?sortDir:0;}
      if(sortK==="dir"){av=dirOf(a);bv=dirOf(b);return av<bv?-sortDir:av>bv?sortDir:0;}
      av=a[sortK];bv=b[sortK]; return (av-bv)*sortDir;
    });
    var html=view.map(function(r){
      var d=dirOf(r);
      return "<tr><td>"+r.gene+"</td><td>"+r.baseMean.toLocaleString("en-GB")+"</td>"+
        "<td>"+r.log2fc.toFixed(2)+"</td><td>"+r.pvalue.toExponential(1)+"</td>"+
        "<td>"+r.padj.toExponential(1)+"</td><td><span class='chip "+d+"'>"+({up:"Up",down:"Down",ns:"n.s."}[d])+"</span></td></tr>";
    }).join("");
    $("tbody").innerHTML=html;
    $("rowcount").textContent="("+view.length+" of "+rows.length+")";
  }
  $("q").addEventListener("input",function(e){filter=e.target.value;render();});
  Array.prototype.forEach.call(document.querySelectorAll("thead th"),function(th){
    th.addEventListener("click",function(){
      var key=th.getAttribute("data-k");
      if(key===sortK){sortDir=-sortDir;} else {sortK=key; sortDir=1;}
      document.querySelectorAll("thead th").forEach(function(t){var b=t.querySelector(".arrow"); if(b)b.remove();});
      var a=document.createElement("span"); a.className="arrow"; a.textContent = sortDir>0?" \u25b2":" \u25bc";
      th.appendChild(a);
      render();
    });
  });

  /* CSV export. Standalone file -> download; embedded viewer (no download rights) -> copy to clipboard. */
  var embedded = (function(){ try{ return window.self !== window.top; }catch(e){ return true; } })();
  function csvText(){
    var head=["gene","gene_id","baseMean","log2FoldChange","pvalue","padj","direction"];
    return [head.join(",")].concat(rows.map(function(r){
      return [r.gene,r.gene_id||"",r.baseMean,r.log2fc,r.pvalue,r.padj,dirOf(r)].join(",");
    })).join("\n");
  }
  if(embedded) $("csv").textContent="Copy CSV";
  $("csv").addEventListener("click",function(){
    var text=csvText(), btn=$("csv"), label=btn.textContent;
    if(embedded){
      var done=function(ok){ btn.textContent=ok?"Copied":"Copy failed"; setTimeout(function(){btn.textContent=label;},1400); };
      try{ navigator.clipboard.writeText(text).then(function(){done(true);},function(){done(false);}); }
      catch(e){ done(false); }
      return;
    }
    try{
      var blob=new Blob([text],{type:"text/csv"});
      var url=URL.createObjectURL(blob);
      var a=document.createElement("a"); a.href=url; a.download="de_results.csv";
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function(){URL.revokeObjectURL(url);},1000);
    }catch(e){}
  });

  /* ---- methods panel: rows from the data when present, plus software and generation time ---- */
  function renderMethods(items){
    var P=D.provenance, el=$("methods"); if(!el) return;
    items=items.slice();
    if(P){
      var v=P.versions||{};
      var pk=Object.keys(v).filter(function(key){return key!=="python";}).map(function(key){return key+" "+v[key];}).join(", ");
      items.push(["Software", pk + (v.python?"; Python "+v.python:"")]);
      items.push(["Generated", (P.generated||"").replace("T"," ")+(P.tool?"  by "+P.tool:"")]);
    }
    if(!items.length) return;
    $("methods-body").innerHTML=items.map(function(r){
      return '<div><div class="m-k">'+r[0]+'</div><div class="m-v">'+r[1]+'</div></div>';}).join("");
    el.hidden=false;
  }
  function legacyMethods(){
    var pr=(D.provenance&&D.provenance.params)||{}, out=[];
    out.push(["Design", D.meta.design||""]);
    out.push(["Counts", (pr.counts_matrix ? "gene counts read from the pipeline count matrix" : "Salmon transcript NumReads summed per gene via the GTF, rounded")+"; genes with < "+(pr.min_count||10)+" total counts dropped"]);
    out.push(["Differential expression", "pyDESeq2, ~"+(pr.condition_col||"condition")+", Wald test, "+(k.shrunk?"apeGLM-style LFC shrinkage, ":"")+"Benjamini-Hochberg; padj < "+PADJ]);
    out.push(["Ordination", "PCA on log2(normalised counts + 1), "+(pr.n_pca_genes||500)+" most variable genes"]);
    out.push(["Correlation", "Spearman on log2 normalised counts, all retained genes"]);
    out.push(["Heatmap", (pr.n_heatmap||28)+" genes with lowest padj, z-scored per gene"]);
    return out;
  }
  function renderHighlights(list){
    if(!list.length) return;
    var box=document.createElement("div"); box.className="highlights";
    box.innerHTML=list.map(function(t){return '<span class="hl">'+t+'</span>';}).join("");
    $("lead").insertAdjacentElement("afterend",box);
  }

  (function(){
    renderMethods(D.methods || legacyMethods());
    var hl=[], S=D.summary;
    if(k.de_genes && S && S.top){
      hl.push("Most significant <b>"+S.top.gene+"</b> padj "+S.top.padj.toExponential(1));
    }
    if(S && S.strongest){
      hl.push("Largest change <b>"+S.strongest.gene+"</b> log2FC "+signed(S.strongest.log2fc,2)+
        " (mean &ge; "+S.strongest.min_basemean+")");
    }
    if(k.de_genes!=null) hl.push("<b>"+fmt(k.de_genes)+"</b> DE genes, "+fmt(k.up)+" up / "+fmt(k.down)+" down");
    if(D.pca&&D.pca.var) hl.push("PC1 explains <b>"+D.pca.var[0]+"%</b> of variance");
    if(S && S.pi0!=null) hl.push("&asymp; <b>"+fmt(S.n_nonnull)+"</b> genes estimated to respond");
    renderHighlights(hl);
  })();

  function drawAll(){ drawPCA(); drawCorr(); drawVolcano(); drawMA(); drawHeat(); drawGenes(); drawPvals(); }
  render(); drawAll();

  function syncBtn(){ $("themeBtn").textContent = currentTheme()==="dark" ? "Light" : "Dark"; }
  syncBtn();
  $("themeBtn").addEventListener("click",function(){
    var next = currentTheme()==="dark" ? "light":"dark";
    document.documentElement.setAttribute("data-theme",next);
    syncBtn(); drawAll();
  });
  if(window.matchMedia){
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change",function(){
      if(!document.documentElement.getAttribute("data-theme")){ syncBtn(); drawAll(); }
    });
  }
})();
