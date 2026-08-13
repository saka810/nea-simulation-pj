"""読み込んだモデルを 3D 表示するビューア（自己完結 HTML + WebGL）。

依存ライブラリなし（標準ライブラリと numpy のみ）。HTML を書き出してブラウザで開く。
matplotlib / pyvista を入れずに済ませているのは、この環境に入っていないため
（`pip install` を前提にしたくない）。将来 GUI を Web ベースにするならこれが土台になる。

表示内容:
  ・三角形要素（辺を描くので分割が見える）
  ・法線ベクトル（矢印）
  ・**法線の裏側を別色で塗る** ← 向きの誤りが一目で分かる
  ・レイヤ別の色分け・表示切り替え
  ・音源 / 受音点
操作: ドラッグで回転、ホイールで拡大縮小、右ドラッグ（または Shift+ドラッグ）で平行移動

使い方:
    python view_model.py ..\\test2.dxf
    python view_model.py ..\\test.dxf --absorption ..\\absorption.csv --out model.html
"""

import argparse
import json
import os
import webbrowser

import numpy as np

import read_dxffile as rd

# レイヤの色（順に割り当て）。彩度を抑えて法線の裏色（赤）と混ざらないようにしている
LAYER_PALETTE = [
    "#5B8FF9", "#61DDAA", "#F6BD16", "#7262FD", "#78D3F8",
    "#9661BC", "#F6903D", "#008685", "#F08BB4", "#B4C6E7",
]


def _hex_to_rgb(code):
    code = code.lstrip("#")
    return [int(code[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def build_payload(model, normal_ratio=0.06):
    """DxfModel から HTML に埋め込む JSON 用の dict を作る。"""
    mesh = model.mesh
    if not mesh:
        raise ValueError("表示できる三角形がありません")

    layers = sorted({t.material for t in mesh})
    layer_color = {name: LAYER_PALETTE[i % len(LAYER_PALETTE)]
                   for i, name in enumerate(layers)}

    lo = model.extents[0]
    hi = model.extents[1]
    diag = float(np.linalg.norm(hi - lo)) or 1.0
    arrow_len = diag * normal_ratio

    positions, normals, colors = [], [], []
    edges = []
    arrows = []
    ranges = []

    # レイヤごとに固めて並べる。こうすると「レイヤの表示切り替え」が
    # 描画範囲を飛ばすだけで済み、シェーダ側に細工が要らない
    for name in layers:
        rgb = _hex_to_rgb(layer_color[name])
        face_start = len(positions) // 3
        edge_start = len(edges) // 3
        arrow_start = len(arrows) // 3

        for t in (m for m in mesh if m.material == name):
            v = np.asarray(t.vertexes, dtype=float)
            n = np.asarray(t.normal, dtype=float)

            for k in range(3):
                positions.extend(v[k].tolist())
                normals.extend(n.tolist())
                colors.extend(rgb)

            # 三角形の辺（分割が見えるように）
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edges.extend(v[a].tolist())
                edges.extend(v[b].tolist())

            # 法線の矢印（重心から法線方向へ）。矢じりは 2 本の短い線で作る
            centre = v.mean(axis=0)
            tip = centre + n * arrow_len
            arrows.extend(centre.tolist())
            arrows.extend(tip.tolist())

            ref = np.array([0.0, 0.0, 1.0])
            if abs(float(np.dot(ref, n))) > 0.9:
                ref = np.array([1.0, 0.0, 0.0])
            side = np.cross(n, ref)
            side = side / (np.linalg.norm(side) or 1.0)
            for sign in (1.0, -1.0):
                back = tip - n * arrow_len * 0.3 + side * arrow_len * 0.15 * sign
                arrows.extend(tip.tolist())
                arrows.extend(back.tolist())

        ranges.append({
            "faceStart": face_start, "faceCount": len(positions) // 3 - face_start,
            "edgeStart": edge_start, "edgeCount": len(edges) // 3 - edge_start,
            "arrowStart": arrow_start, "arrowCount": len(arrows) // 3 - arrow_start,
        })

    return {
        "positions": positions,
        "normals": normals,
        "colors": colors,
        "edges": edges,
        "arrows": arrows,
        "layers": [dict({"name": name, "color": layer_color[name],
                         "count": model.layer_counts.get(name, 0)}, **ranges[i])
                   for i, name in enumerate(layers)],
        "source": [p.tolist() for p in model.source_points],
        "receiver": [p.tolist() for p in model.receiver_points],
        "bbox": [lo.tolist(), hi.tolist()],
        "arrowLength": arrow_len,
        "summary": model.summary(),
        "triangleCount": len(mesh),
    }


HTML_TEMPLATE = r"""<title>__TITLE__</title>
<style>
:root{
  --bg:#0E1319; --panel:#161D26; --panel-2:#1D2732; --rule:#2A3542;
  --ink:#E4EAF1; --ink-2:#A8B5C4; --muted:#7C8B9C; --accent:#F6BD16;
  --front:#5B8FF9; --back:#C2453C;
  --mono:ui-monospace,"SFMono-Regular",Consolas,monospace;
  --sans:"Hiragino Kaku Gothic ProN","Yu Gothic","Meiryo",system-ui,sans-serif;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);overflow:hidden}
#wrap{position:fixed;inset:0;display:flex}
canvas{flex:1;display:block;cursor:grab;background:var(--bg)}
canvas.dragging{cursor:grabbing}

#panel{
  width:300px;flex:none;background:var(--panel);border-left:1px solid var(--rule);
  overflow-y:auto;padding:18px 18px 28px;font-size:13px;line-height:1.7;
}
#panel h1{font-size:15px;margin:0 0 4px;font-weight:700}
#panel .sub{color:var(--muted);font-size:11.5px;margin:0 0 18px;word-break:break-all}
#panel h2{
  font-size:10.5px;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);
  margin:22px 0 9px;padding-bottom:6px;border-bottom:1px solid var(--rule);font-weight:700;
}
#panel h2:first-of-type{margin-top:0}
label.row{display:flex;align-items:center;gap:9px;padding:3px 0;cursor:pointer;user-select:none}
label.row:hover{color:#fff}
input[type=checkbox]{accent-color:var(--accent);width:14px;height:14px;flex:none;margin:0}
.swatch{width:11px;height:11px;border-radius:2px;flex:none;border:1px solid rgba(255,255,255,.25)}
.count{margin-left:auto;color:var(--muted);font-family:var(--mono);font-size:11px}
.legend{display:flex;align-items:center;gap:9px;padding:3px 0}
.btns{display:grid;grid-template-columns:1fr 1fr;gap:6px;margin-top:4px}
button{
  font-family:var(--sans);font-size:12px;padding:7px 4px;border-radius:3px;cursor:pointer;
  background:var(--panel-2);color:var(--ink);border:1px solid var(--rule);
}
button:hover{background:#25313E;border-color:#3A4A5C}
button:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
pre{
  margin:0;font-family:var(--mono);font-size:10.5px;line-height:1.65;color:var(--ink-2);
  white-space:pre-wrap;word-break:break-word;
}
.hint{color:var(--muted);font-size:11px;line-height:1.6}
kbd{
  font-family:var(--mono);font-size:10px;background:var(--panel-2);border:1px solid var(--rule);
  border-radius:2px;padding:1px 4px;color:var(--ink);
}
#readout{
  position:fixed;left:14px;bottom:12px;font-family:var(--mono);font-size:11px;
  color:var(--muted);pointer-events:none;
}
</style>

<div id="wrap">
  <canvas id="gl"></canvas>
  <div id="panel">
    <h1>モデルビューア</h1>
    <p class="sub" id="fname"></p>

    <h2>表示</h2>
    <label class="row"><input type="checkbox" id="cbFaces" checked><span>面</span></label>
    <label class="row"><input type="checkbox" id="cbEdges" checked><span>三角形の辺</span></label>
    <label class="row"><input type="checkbox" id="cbNormals" checked><span>法線ベクトル</span></label>
    <label class="row"><input type="checkbox" id="cbBack" checked><span>法線の裏側を赤で塗る</span></label>
    <label class="row"><input type="checkbox" id="cbPoints" checked><span>音源・受音点</span></label>
    <label class="row"><input type="checkbox" id="cbAxes" checked><span>座標軸</span></label>

    <h2>法線の向き</h2>
    <div class="legend"><span class="swatch" style="background:var(--front)"></span>
      <span>表（法線がこちらを向く）</span></div>
    <div class="legend"><span class="swatch" style="background:var(--back)"></span>
      <span>裏（音線が通り抜ける側）</span></div>
    <p class="hint">面は片側だけ反射します。赤く見える側から来た音線はその面を通り抜けます。</p>

    <h2>レイヤ</h2>
    <div id="layers"></div>

    <h2>視点</h2>
    <div class="btns">
      <button data-view="iso">等角</button>
      <button data-view="top">上（平面）</button>
      <button data-view="front">正面 (+Y)</button>
      <button data-view="side">側面 (+X)</button>
    </div>
    <p class="hint" style="margin-top:10px">
      <kbd>ドラッグ</kbd> 回転　<kbd>ホイール</kbd> 拡大縮小<br>
      <kbd>右ドラッグ</kbd> または <kbd>Shift</kbd>+ドラッグ 平行移動
    </p>

    <h2>読み込み結果</h2>
    <pre id="summary"></pre>
  </div>
</div>
<div id="readout"></div>

<script>
const DATA = __DATA__;

// ---------- 行列（列優先。WebGL の並びに合わせる） ----------
function mIdent(){return [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];}
function mMul(a,b){
  const o=new Array(16);
  for(let c=0;c<4;c++)for(let r=0;r<4;r++){
    let s=0; for(let k=0;k<4;k++) s+=a[k*4+r]*b[c*4+k];
    o[c*4+r]=s;
  }
  return o;
}
function mPerspective(fovy,aspect,near,far){
  const f=1/Math.tan(fovy/2), nf=1/(near-far);
  return [f/aspect,0,0,0, 0,f,0,0, 0,0,(far+near)*nf,-1, 0,0,2*far*near*nf,0];
}
function vSub(a,b){return [a[0]-b[0],a[1]-b[1],a[2]-b[2]];}
function vCross(a,b){return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];}
function vNorm(a){const l=Math.hypot(a[0],a[1],a[2])||1;return [a[0]/l,a[1]/l,a[2]/l];}
function vDot(a,b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
function mLookAt(eye,center,up){
  const f=vNorm(vSub(center,eye)), s=vNorm(vCross(f,up)), u=vCross(s,f);
  return [s[0],u[0],-f[0],0, s[1],u[1],-f[1],0, s[2],u[2],-f[2],0,
          -vDot(s,eye),-vDot(u,eye),vDot(f,eye),1];
}

// ---------- WebGL の下準備 ----------
const canvas=document.getElementById("gl");
const gl=canvas.getContext("webgl",{antialias:true,alpha:false});
if(!gl){document.getElementById("summary").textContent="WebGL が使えません。";}

function shader(type,src){
  const s=gl.createShader(type); gl.shaderSource(s,src); gl.compileShader(s);
  if(!gl.getShaderParameter(s,gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}
function program(vs,fs){
  const p=gl.createProgram();
  gl.attachShader(p,shader(gl.VERTEX_SHADER,vs));
  gl.attachShader(p,shader(gl.FRAGMENT_SHADER,fs));
  gl.linkProgram(p);
  if(!gl.getProgramParameter(p,gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}

// 面：法線が視線側を向いていれば「表」、逆なら「裏」として色を変える。
// gl_FrontFacing（頂点の巻き順）ではなく **保持している法線** で判定するのが要点
// （読み込み側は巻き順を変えずに法線だけ反転することがあるため）。
const progFace=program(`
attribute vec3 aPos; attribute vec3 aNormal; attribute vec3 aColor;
uniform mat4 uMVP; uniform vec3 uEye;
varying vec3 vN; varying vec3 vC; varying vec3 vToEye;
void main(){
  gl_Position = uMVP * vec4(aPos,1.0);
  vN=aNormal; vC=aColor; vToEye=uEye-aPos;
}`,`
precision mediump float;
varying vec3 vN; varying vec3 vC; varying vec3 vToEye;
uniform float uBackMark;
void main(){
  vec3 n=normalize(vN); vec3 e=normalize(vToEye);
  float facing=dot(n,e);
  float lambert=abs(facing)*0.72+0.28;
  vec3 base=vC;
  if(facing<0.0 && uBackMark>0.5){ base=vec3(0.76,0.27,0.24); }
  gl_FragColor=vec4(base*lambert,1.0);
}`);

const progLine=program(`
attribute vec3 aPos; uniform mat4 uMVP;
void main(){ gl_Position=uMVP*vec4(aPos,1.0); }`,`
precision mediump float; uniform vec4 uColor;
void main(){ gl_FragColor=uColor; }`);

const progPoint=program(`
attribute vec3 aPos; uniform mat4 uMVP; uniform float uSize;
void main(){ gl_Position=uMVP*vec4(aPos,1.0); gl_PointSize=uSize; }`,`
precision mediump float; uniform vec4 uColor;
void main(){
  vec2 d=gl_PointCoord-vec2(0.5);
  if(dot(d,d)>0.25) discard;
  gl_FragColor=uColor;
}`);

function buffer(arr){
  const b=gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER,b);
  gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(arr),gl.STATIC_DRAW);
  return b;
}

const bPos=buffer(DATA.positions), bNrm=buffer(DATA.normals),
      bCol=buffer(DATA.colors),
      bEdge=buffer(DATA.edges), bArrow=buffer(DATA.arrows);

const srcFlat=[].concat(...DATA.source), rcvFlat=[].concat(...DATA.receiver);
const bSrc=srcFlat.length?buffer(srcFlat):null;
const bRcv=rcvFlat.length?buffer(rcvFlat):null;

// 座標軸（原点からモデルサイズの 15%）
const lo=DATA.bbox[0], hi=DATA.bbox[1];
const centre=[(lo[0]+hi[0])/2,(lo[1]+hi[1])/2,(lo[2]+hi[2])/2];
const diag=Math.hypot(hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2])||1;
const axLen=diag*0.15;
const bAxX=buffer([lo[0],lo[1],lo[2], lo[0]+axLen,lo[1],lo[2]]);
const bAxY=buffer([lo[0],lo[1],lo[2], lo[0],lo[1]+axLen,lo[2]]);
const bAxZ=buffer([lo[0],lo[1],lo[2], lo[0],lo[1],lo[2]+axLen]);

// ---------- カメラ ----------
const cam={az:-Math.PI*0.32, el:Math.PI*0.22, dist:diag*2.0,
           target:centre.slice()};

function eyePos(){
  const ce=Math.cos(cam.el), se=Math.sin(cam.el);
  return [cam.target[0]+cam.dist*ce*Math.cos(cam.az),
          cam.target[1]+cam.dist*ce*Math.sin(cam.az),
          cam.target[2]+cam.dist*se];
}

function setView(kind){
  cam.target=centre.slice(); cam.dist=diag*2.0;
  if(kind==="iso"){ cam.az=-Math.PI*0.32; cam.el=Math.PI*0.22; }
  // 真上（el=90°）だと up ベクトル [0,0,1] と視線が平行になり lookAt が退化するので
  // わずかに傾ける（見た目はほぼ真上のまま）
  if(kind==="top"){ cam.az=-Math.PI/2; cam.el=Math.PI/2-0.02; }
  if(kind==="front"){ cam.az=-Math.PI/2; cam.el=0.0; }
  if(kind==="side"){ cam.az=0.0; cam.el=0.0; }
  draw();
}

// ---------- 描画 ----------
// レイヤの表示切り替えは、レイヤ順に並べた頂点の描画範囲を飛ばすだけで実現する
const visible=DATA.layers.map(()=>true);

function bindAttr(prog,name,buf,size){
  const loc=gl.getAttribLocation(prog,name);
  if(loc<0) return;
  gl.bindBuffer(gl.ARRAY_BUFFER,buf);
  gl.enableVertexAttribArray(loc);
  gl.vertexAttribPointer(loc,size,gl.FLOAT,false,0,0);
}

function drawLines(buf,start,count,color,mvp){
  if(count<=0) return;
  gl.useProgram(progLine);
  gl.uniformMatrix4fv(gl.getUniformLocation(progLine,"uMVP"),false,mvp);
  gl.uniform4fv(gl.getUniformLocation(progLine,"uColor"),color);
  bindAttr(progLine,"aPos",buf,3);
  gl.drawArrays(gl.LINES,start,count);
}

function drawPoints(buf,count,color,size,mvp){
  gl.useProgram(progPoint);
  gl.uniformMatrix4fv(gl.getUniformLocation(progPoint,"uMVP"),false,mvp);
  gl.uniform4fv(gl.getUniformLocation(progPoint,"uColor"),color);
  gl.uniform1f(gl.getUniformLocation(progPoint,"uSize"),size);
  bindAttr(progPoint,"aPos",buf,3);
  gl.drawArrays(gl.POINTS,0,count);
}

function draw(){
  const dpr=Math.min(window.devicePixelRatio||1,2);
  const w=canvas.clientWidth, h=canvas.clientHeight;
  canvas.width=w*dpr; canvas.height=h*dpr;
  gl.viewport(0,0,canvas.width,canvas.height);
  gl.clearColor(0.055,0.075,0.098,1);
  gl.enable(gl.DEPTH_TEST);
  gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);

  const eye=eyePos();
  const proj=mPerspective(Math.PI/4, w/h, diag*0.002, diag*20);
  const view=mLookAt(eye,cam.target,[0,0,1]);
  const mvp=mMul(proj,view);

  if(document.getElementById("cbFaces").checked){
    gl.useProgram(progFace);
    gl.uniformMatrix4fv(gl.getUniformLocation(progFace,"uMVP"),false,mvp);
    gl.uniform3fv(gl.getUniformLocation(progFace,"uEye"),eye);
    gl.uniform1f(gl.getUniformLocation(progFace,"uBackMark"),
                 document.getElementById("cbBack").checked?1:0);
    bindAttr(progFace,"aPos",bPos,3);
    bindAttr(progFace,"aNormal",bNrm,3);
    bindAttr(progFace,"aColor",bCol,3);
    gl.enable(gl.POLYGON_OFFSET_FILL);
    gl.polygonOffset(1.0,1.0);
    DATA.layers.forEach((L,i)=>{
      if(visible[i] && L.faceCount>0) gl.drawArrays(gl.TRIANGLES,L.faceStart,L.faceCount);
    });
    gl.disable(gl.POLYGON_OFFSET_FILL);
  }
  if(document.getElementById("cbEdges").checked)
    DATA.layers.forEach((L,i)=>{
      if(visible[i]) drawLines(bEdge,L.edgeStart,L.edgeCount,[0.86,0.90,0.95,1.0],mvp);
    });
  if(document.getElementById("cbNormals").checked)
    DATA.layers.forEach((L,i)=>{
      if(visible[i]) drawLines(bArrow,L.arrowStart,L.arrowCount,[0.965,0.741,0.086,1.0],mvp);
    });
  if(document.getElementById("cbAxes").checked){
    drawLines(bAxX,0,2,[0.90,0.35,0.32,1],mvp);
    drawLines(bAxY,0,2,[0.38,0.80,0.45,1],mvp);
    drawLines(bAxZ,0,2,[0.40,0.60,0.95,1],mvp);
  }
  if(document.getElementById("cbPoints").checked){
    if(bSrc) drawPoints(bSrc,DATA.source.length,[0.95,0.30,0.25,1],14*dpr,mvp);
    if(bRcv) drawPoints(bRcv,DATA.receiver.length,[0.35,0.72,0.98,1],14*dpr,mvp);
  }

  document.getElementById("readout").textContent=
    "方位 "+(cam.az*180/Math.PI).toFixed(0)+"°  仰角 "+(cam.el*180/Math.PI).toFixed(0)+
    "°  距離 "+cam.dist.toFixed(2)+" m";
}

// ---------- 操作 ----------
let drag=null;
canvas.addEventListener("mousedown",e=>{
  drag={x:e.clientX,y:e.clientY,pan:(e.button===2||e.shiftKey)};
  canvas.classList.add("dragging"); e.preventDefault();
});
window.addEventListener("mouseup",()=>{drag=null;canvas.classList.remove("dragging");});
window.addEventListener("mousemove",e=>{
  if(!drag) return;
  const dx=e.clientX-drag.x, dy=e.clientY-drag.y;
  drag.x=e.clientX; drag.y=e.clientY;
  if(drag.pan){
    const eye=eyePos();
    const fwd=vNorm(vSub(cam.target,eye));
    const right=vNorm(vCross(fwd,[0,0,1]));
    const up=vCross(right,fwd);
    const k=cam.dist*0.0016;
    for(let i=0;i<3;i++) cam.target[i]+= -right[i]*dx*k + up[i]*dy*k;
  }else{
    cam.az-=dx*0.008;
    cam.el=Math.max(-Math.PI/2+0.01,Math.min(Math.PI/2-0.01,cam.el+dy*0.008));
  }
  draw();
});
canvas.addEventListener("contextmenu",e=>e.preventDefault());
canvas.addEventListener("wheel",e=>{
  cam.dist=Math.max(diag*0.05,Math.min(diag*15,cam.dist*(1+Math.sign(e.deltaY)*0.12)));
  draw(); e.preventDefault();
},{passive:false});
window.addEventListener("resize",draw);

// ---------- UI ----------
document.getElementById("fname").textContent=DATA.title||"";
document.getElementById("summary").textContent=DATA.summary;
const layerBox=document.getElementById("layers");
DATA.layers.forEach((L,i)=>{
  const lab=document.createElement("label");
  lab.className="row";
  lab.innerHTML='<input type="checkbox" checked data-layer="'+i+'">'
    +'<span class="swatch" style="background:'+L.color+'"></span>'
    +'<span>'+L.name+'</span><span class="count">'+L.count+'</span>';
  layerBox.appendChild(lab);
});
layerBox.addEventListener("change",e=>{
  const i=parseInt(e.target.dataset.layer,10);
  if(!isNaN(i)) visible[i]=e.target.checked;
  draw();
});
document.querySelectorAll("#panel input[type=checkbox]:not([data-layer])")
  .forEach(cb=>cb.addEventListener("change",draw));
document.querySelectorAll("button[data-view]")
  .forEach(b=>b.addEventListener("click",()=>setView(b.dataset.view)));

// レイアウト確定後に描く（canvas.clientWidth が 0 のまま描かないように）
requestAnimationFrame(draw);
</script>
"""


def export_html(model, out_path, title="", subtitle=""):
    """DxfModel を自己完結 HTML に書き出す。

    title    … ブラウザのタブに出る名前
    subtitle … パネル上部に出す説明（元の DXF ファイル名など）
    """
    payload = build_payload(model)
    payload["title"] = subtitle or title
    html = (HTML_TEMPLATE
            .replace("__DATA__", json.dumps(payload, ensure_ascii=False))
            .replace("__TITLE__", title or "モデルビューア"))
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    return out_path


def view(dxf_path, out_path=None, absorption=None, unit=None,
         orient_normals="cad", open_browser=True):
    """DXF を読み込んで 3D ビューアの HTML を書き出し、ブラウザで開く。"""
    model = rd.read_model(dxf_path, unit=unit, absorption_table=absorption,
                          orient_normals=orient_normals)
    base = os.path.splitext(os.path.basename(dxf_path))[0]
    if out_path is None:
        out_path = os.path.join(os.path.dirname(os.path.abspath(dxf_path)),
                                base + "_view.html")
    export_html(model, out_path, title=f"{base} モデルビューア",
                subtitle=os.path.basename(dxf_path))
    print(f"\n[view_model] 書き出しました: {out_path}")
    if open_browser:
        webbrowser.open("file:///" + os.path.abspath(out_path).replace("\\", "/"))
    return out_path, model


def main():
    p = argparse.ArgumentParser(description="DXF モデルの 3D ビューア（HTML を書き出す）")
    p.add_argument("dxf", help="DXF ファイル")
    p.add_argument("--out", help="出力する HTML のパス")
    p.add_argument("--absorption", help="吸音率 CSV")
    p.add_argument("--unit", help="'mm' / 'm' など。省略すると $INSUNITS から自動判定")
    p.add_argument("--orient-normals", default="cad", choices=["cad", "flip", "shells"])
    p.add_argument("--no-open", action="store_true", help="ブラウザを開かない")
    a = p.parse_args()
    view(a.dxf, out_path=a.out, absorption=a.absorption, unit=a.unit,
         orient_normals=a.orient_normals, open_browser=not a.no_open)


if __name__ == "__main__":
    main()
