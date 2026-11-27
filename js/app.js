// Interface : dessin main levée, outils, affichage des résultats
(function () {
  "use strict";
  const FT = window.FT;
  const $ = (id) => document.getElementById(id);
  const KEY = "fieldtool.project.v1";

  const defaults = () => ({ edges: [], labels: [], customMats: [], grid: 5, gridReal: 10, drawMode: "grid", depth: 20, physics: "mag", airMargin: 1, meshN: 40, bcElec: "V0" });
  let S = defaults();
  let view = { x: 30, y: 30, z: 6 }; // centre (mm) et zoom (px/mm)
  let tool = "draw", sel = -1, mesh = null, result = null, post = null, probe = null, bg = null;
  let armed = null; // matériau choisi dans la palette (pose directe avec l'outil Matériau)
  let act = null, penSeen = false, meshTimer = 0, meshError = "";
  const disp = { qty: "B", lines: true, arrows: false, mesh: false, dims: true, max: null };
  const undo = [], redo = [];

  const wrap = $("wrap"), glc = $("gl"), cv = $("cv"), ctx = cv.getContext("2d");
  let W = 100, H = 100, dpr = 1;

  // ---------- utilitaires ----------
  const w2s = (x, y) => [W / 2 + (x - view.x) * view.z, H / 2 - (y - view.y) * view.z];
  const s2w = (sx, sy) => ({ x: view.x + (sx - W / 2) / view.z, y: view.y - (sy - H / 2) / view.z });
  const snap = (p) => ({ x: Math.round(p.x / S.grid) * S.grid, y: Math.round(p.y / S.grid) * S.grid });
  // accrochage selon le mode : grille, ou sommets/arêtes voisins en mode patate
  const snapPt = (p, skip) => (S.drawMode === "free" ? FT.snapFree(p, S.edges, 12 / view.z, skip) : snap(p));
  const num = (v, d) => { const n = parseFloat(String(v).replace(",", ".")); return isFinite(n) ? n : d; };
  function fmt(v, unit) {
    if (!isFinite(v)) return "–";
    const a = Math.abs(v);
    if (a === 0) return "0 " + unit;
    const pre = [[1e9, "G"], [1e6, "M"], [1e3, "k"], [1, ""], [1e-3, "m"], [1e-6, "µ"], [1e-9, "n"], [1e-12, "p"], [1e-15, "f"]];
    for (const [f, p] of pre) if (a >= f * 0.9995) return (v / f).toPrecision(4) + " " + p + unit;
    return v.toExponential(3) + " " + unit;
  }
  const mm = (v) => +v.toFixed(3) + " mm";
  let toastTimer = 0;
  function toast(msg, ms) {
    const t = $("toast");
    t.textContent = msg; t.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => (t.hidden = true), ms || 2600);
  }
  // Dialogue intégré à la page (confirm/prompt du navigateur ne sont pas toujours disponibles).
  // o : {msg, input?, text?, readOnly?, buttons:[[libellé, valeur | fonction(done), principal?]]}
  function ask(o) {
    return new Promise((res) => {
      const d = $("ask"), i = $("askIn"), t = $("askTxt"), bx = $("askBtns");
      $("askMsg").textContent = o.msg;
      i.hidden = o.input == null; if (o.input != null) i.value = o.input;
      t.hidden = o.text == null; if (o.text != null) t.value = o.text;
      t.readOnly = !!o.readOnly;
      const done = (v) => { d.close(); res(v); };
      d.oncancel = () => res(null);
      bx.innerHTML = ""; bx.className = o.list ? "btns list" : "btns end";
      for (const [label, val, primary] of o.buttons) {
        const b = document.createElement("button");
        b.type = "button"; b.textContent = label;
        if (primary) b.className = "primary";
        b.onclick = () => (typeof val === "function" ? val(done, b) : done(val === true ? (o.input != null ? i.value : o.text != null ? t.value : true) : val));
        bx.appendChild(b);
      }
      d.showModal();
    });
  }
  const hexRgb = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  const STOPS = [[48, 18, 59], [70, 107, 227], [40, 187, 236], [50, 242, 152], [164, 252, 60], [238, 207, 58], [251, 126, 33], [209, 48, 5], [122, 4, 3]];
  function cmap(t, out, o) {
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    const u = t * 8, i = Math.min(7, Math.floor(u)), f = u - i, a = STOPS[i], b = STOPS[i + 1];
    out[o] = a[0] + (b[0] - a[0]) * f; out[o + 1] = a[1] + (b[1] - a[1]) * f; out[o + 2] = a[2] + (b[2] - a[2]) * f; out[o + 3] = 255;
  }

  // ---------- WebGL (remplissage des triangles) ----------
  let gl = null, glProg = null, glPos = null, glCol = null, glN = 0, glT = null, colArr = null;
  function initGL() {
    try { gl = glc.getContext("webgl", { antialias: true, alpha: false }); } catch (e) { gl = null; }
    if (!gl) return;
    const sh = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); return s; };
    glProg = gl.createProgram();
    gl.attachShader(glProg, sh(gl.VERTEX_SHADER, "attribute vec2 p;attribute vec4 c;uniform vec4 t;varying vec4 v;void main(){v=c;gl_Position=vec4(p.x*t.x+t.z,p.y*t.y+t.w,0.,1.);}"));
    gl.attachShader(glProg, sh(gl.FRAGMENT_SHADER, "precision mediump float;varying vec4 v;void main(){gl_FragColor=v;}"));
    gl.linkProgram(glProg);
    if (!gl.getProgramParameter(glProg, gl.LINK_STATUS)) { gl = null; return; }
    glPos = gl.createBuffer(); glCol = gl.createBuffer();
    glT = gl.getUniformLocation(glProg, "t");
  }
  function uploadMesh() {
    glN = 0;
    if (!mesh) return;
    const nt = mesh.tri.length / 3;
    const pos = new Float32Array(6 * nt);
    for (let k = 0; k < 3 * nt; k++) { pos[2 * k] = mesh.x[mesh.tri[k]]; pos[2 * k + 1] = mesh.y[mesh.tri[k]]; }
    colArr = new Uint8Array(12 * nt);
    glN = 3 * nt;
    if (gl) { gl.bindBuffer(gl.ARRAY_BUFFER, glPos); gl.bufferData(gl.ARRAY_BUFFER, pos, gl.STATIC_DRAW); }
  }
  function uploadColors() {
    if (!mesh) return;
    const nt = mesh.tri.length / 3, c = colArr;
    if (result && post) {
      const v = post.val, lo = post.lo, sp = post.hi - post.lo || 1;
      for (let k = 0; k < 3 * nt; k++) cmap((v[k] - lo) / sp, c, 4 * k);
    } else {
      const rc = [];
      for (let r = 0; r < mesh.nReg; r++) { const l = mesh.regLabel[r]; rc.push(l ? hexRgb(FT.getMaterial(S, l.mat).color) : [255, 255, 255]); }
      for (let e = 0; e < nt; e++) {
        const q = rc[mesh.reg[e]];
        for (let a = 0; a < 3; a++) { const o = 12 * e + 4 * a; c[o] = q[0]; c[o + 1] = q[1]; c[o + 2] = q[2]; c[o + 3] = 255; }
      }
    }
    if (gl) { gl.bindBuffer(gl.ARRAY_BUFFER, glCol); gl.bufferData(gl.ARRAY_BUFFER, c, gl.DYNAMIC_DRAW); }
  }
  function drawGL() {
    if (!gl) return;
    gl.viewport(0, 0, glc.width, glc.height);
    gl.clearColor(0.933, 0.945, 0.957, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    if (!mesh || !glN) return;
    gl.useProgram(glProg);
    const sx = (2 * view.z) / W, sy = (2 * view.z) / H;
    gl.uniform4f(glT, sx, sy, -view.x * sx, -view.y * sy);
    const lp = gl.getAttribLocation(glProg, "p"), lc = gl.getAttribLocation(glProg, "c");
    gl.bindBuffer(gl.ARRAY_BUFFER, glPos); gl.enableVertexAttribArray(lp); gl.vertexAttribPointer(lp, 2, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, glCol); gl.enableVertexAttribArray(lc); gl.vertexAttribPointer(lc, 4, gl.UNSIGNED_BYTE, true, 0, 0);
    gl.drawArrays(gl.TRIANGLES, 0, glN);
  }

  // ---------- post-traitement pour l'affichage ----------
  const QTY = {
    mag: [["B", "|B| – induction (T)"], ["H", "|H| – champ (A/m)"]],
    elec: [["V", "V – potentiel (V)"], ["E", "|E| – champ (V/mm)"]],
  };
  function buildPost() {
    post = null;
    if (!result || !mesh) return;
    const nt = mesh.tri.length / 3, r = result;
    let val, lo = 0, hi, unit, title, elemVals;
    if (r.type === "mag") {
      elemVals = disp.qty === "H" ? r.Hm : r.Bm;
      unit = disp.qty === "H" ? "A/m" : "T"; title = disp.qty === "H" ? "|H|" : "|B|";
      val = FT.smoothPerRegion(mesh, elemVals, r.area);
    } else if (disp.qty === "E") {
      elemVals = Float64Array.from(r.Em, (v) => v * 1e-3);
      unit = "V/mm"; title = "|E|";
      val = FT.smoothPerRegion(mesh, elemVals, r.area);
    } else {
      val = new Float32Array(3 * nt);
      lo = Infinity; hi = -Infinity;
      for (let k = 0; k < 3 * nt; k++) { const v = r.V[mesh.tri[k]]; val[k] = v; if (v < lo) lo = v; if (v > hi) hi = v; }
      unit = "V"; title = "V";
    }
    if (elemVals) {
      const s = Float64Array.from(elemVals).sort();
      hi = s[Math.min(s.length - 1, Math.floor(s.length * 0.995))] || 1;
      if (disp.max > 0) hi = disp.max;
    }
    if (!(hi > lo)) hi = lo + 1;
    const vx = r.type === "mag" ? r.Bx : r.Ex, vy = r.type === "mag" ? r.By : r.Ey;
    post = {
      val, lo, hi, unit, title,
      sx: FT.smoothPerRegion(mesh, vx, r.area), sy: FT.smoothPerRegion(mesh, vy, r.area),
      lines: FT.contours(mesh, r.type === "mag" ? r.A : r.V, 26),
      stats: {},
    };
  }

  // ---------- dessin ----------
  function resize() {
    const r = wrap.getBoundingClientRect();
    W = Math.max(50, r.width); H = Math.max(50, r.height);
    dpr = window.devicePixelRatio || 1;
    for (const c of [glc, cv]) { c.width = Math.round(W * dpr); c.height = Math.round(H * dpr); }
    draw();
  }
  let rafPending = false;
  function draw() {
    if (rafPending) return;
    rafPending = true;
    requestAnimationFrame(() => { rafPending = false; render(); });
  }
  function niceStep(v) {
    const p = Math.pow(10, Math.floor(Math.log10(v))), f = v / p;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
  }
  function render() {
    drawGL();
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    const a = s2w(0, H), b = s2w(W, 0);
    if (!gl && mesh) { // repli sans WebGL
      const nt = mesh.tri.length / 3;
      for (let e = 0; e < nt; e++) {
        ctx.beginPath();
        for (let k = 0; k < 3; k++) { const p = w2s(mesh.x[mesh.tri[3 * e + k]], mesh.y[mesh.tri[3 * e + k]]); k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]); }
        const o = 12 * e;
        ctx.fillStyle = "rgb(" + ((colArr[o] + colArr[o + 4] + colArr[o + 8]) / 3 | 0) + "," + ((colArr[o + 1] + colArr[o + 5] + colArr[o + 9]) / 3 | 0) + "," + ((colArr[o + 2] + colArr[o + 6] + colArr[o + 10]) / 3 | 0) + ")";
        ctx.fill();
      }
    }
    if (bg) {
      const p = w2s(bg.x, bg.y);
      ctx.globalAlpha = bg.alpha;
      ctx.drawImage(bg.img, p[0], p[1], bg.img.width * bg.mm * view.z, bg.img.height * bg.mm * view.z);
      ctx.globalAlpha = 1;
    }
    // grille
    if (!result) {
      // grille réelle : pas fixe en mm, élargi en 1-2-5 si les lignes deviennent trop serrées
      let sr = S.gridReal, kf = 0;
      while (sr * view.z < 24) sr *= [2, 2.5, 2][kf++ % 3];
      const lines = (st, style) => {
        ctx.beginPath(); ctx.strokeStyle = style;
        for (let i = Math.ceil(a.x / st); i * st <= b.x; i++) { const x = Math.round(w2s(i * st, 0)[0]) + 0.5; ctx.moveTo(x, 0); ctx.lineTo(x, H); }
        for (let j = Math.ceil(a.y / st); j * st <= b.y; j++) { const y = Math.round(w2s(0, j * st)[1]) + 0.5; ctx.moveTo(0, y); ctx.lineTo(W, y); }
        ctx.stroke();
      };
      ctx.lineWidth = 1;
      // grille de dessin (accrochage) : masquée quand elle est trop fine à l'écran
      if (S.drawMode !== "free" && S.grid * view.z >= 6 && S.grid < sr) lines(S.grid, "rgba(60,80,105,.11)");
      lines(sr, "rgba(20,99,214,.34)");
      const nf = (v) => String(+v.toPrecision(10));
      const every = Math.max(1, Math.ceil(46 / (sr * view.z)));
      ctx.font = "11px -apple-system, Segoe UI, sans-serif"; ctx.fillStyle = "rgba(20,70,150,.85)";
      ctx.textBaseline = "top"; ctx.textAlign = "left";
      for (let i = Math.ceil(a.x / sr); i * sr <= b.x; i++) if (i % every === 0) ctx.fillText(nf(i * sr), w2s(i * sr, 0)[0] + 3, 3);
      for (let j = Math.ceil(a.y / sr); j * sr <= b.y; j++) { const y = w2s(0, j * sr)[1]; if (j % every === 0 && y > 20) ctx.fillText(nf(j * sr), 3, y + 2); }
      ctx.textBaseline = "bottom"; ctx.textAlign = "right";
      ctx.fillText("Grille réelle : " + nf(sr) + " mm" + (sr !== S.gridReal ? " (réglée à " + nf(S.gridReal) + ")" : "") + "   ·   dessin : " + nf(S.grid) + " mm", W - 8, H - 6);
      ctx.textAlign = "left";
    }
    if (mesh) {
      const p0 = w2s(mesh.box.x0, mesh.box.y1), p1 = w2s(mesh.box.x1, mesh.box.y0);
      ctx.setLineDash([6, 4]); ctx.strokeStyle = "#8894a1"; ctx.lineWidth = 1;
      ctx.strokeRect(p0[0], p0[1], p1[0] - p0[0], p1[1] - p0[1]);
      ctx.setLineDash([]);
      if (disp.mesh) {
        ctx.beginPath();
        const nt = mesh.tri.length / 3;
        for (let e = 0; e < nt; e++) {
          const q0 = w2s(mesh.x[mesh.tri[3 * e]], mesh.y[mesh.tri[3 * e]]), q1 = w2s(mesh.x[mesh.tri[3 * e + 1]], mesh.y[mesh.tri[3 * e + 1]]), q2 = w2s(mesh.x[mesh.tri[3 * e + 2]], mesh.y[mesh.tri[3 * e + 2]]);
          ctx.moveTo(q0[0], q0[1]); ctx.lineTo(q1[0], q1[1]); ctx.lineTo(q2[0], q2[1]); ctx.closePath();
        }
        ctx.strokeStyle = result ? "rgba(255,255,255,.25)" : "rgba(40,60,90,.25)"; ctx.lineWidth = 0.5; ctx.stroke();
      }
    }
    if (result && post) {
      if (disp.lines) {
        const L = post.lines;
        ctx.beginPath();
        for (let i = 0; i < L.length; i += 4) { const q0 = w2s(L[i], L[i + 1]), q1 = w2s(L[i + 2], L[i + 3]); ctx.moveTo(q0[0], q0[1]); ctx.lineTo(q1[0], q1[1]); }
        ctx.strokeStyle = "rgba(255,255,255,.75)"; ctx.lineWidth = 1; ctx.stroke();
      }
      if (disp.arrows) drawArrows();
    }
    // arêtes
    ctx.lineCap = "round"; ctx.lineJoin = "round";
    const erasing = act && act.type === "erase" ? act.hit : null;
    S.edges.forEach((e, i) => {
      const pts = FT.edgePoints(e, 4);
      ctx.beginPath();
      pts.forEach((p, k) => { const q = w2s(p[0], p[1]); k ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); });
      ctx.strokeStyle = erasing && erasing.has(i) ? "#e11d48" : i === sel ? "#1463d6" : result ? "#0b1118" : "#1c2733";
      ctx.lineWidth = i === sel ? 3.5 : 2;
      ctx.stroke();
    });
    if (tool === "edit" || tool === "draw") {
      ctx.fillStyle = "#1c2733";
      const seen = new Set();
      for (const e of S.edges)
        for (const p of [[e.ax, e.ay], [e.bx, e.by]]) {
          const k = p[0] + "," + p[1];
          if (seen.has(k)) continue;
          seen.add(k);
          const q = w2s(p[0], p[1]);
          ctx.beginPath(); ctx.arc(q[0], q[1], tool === "edit" ? 4.5 : 2.5, 0, 6.3); ctx.fill();
        }
    }
    if (disp.dims && !result) {
      ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      for (const e of S.edges) {
        const len = FT.edgeLength(e);
        if (len * view.z < 46) continue;
        let mx, my, txt;
        if (e.type === "arc") { const ap = FT.arcParams(e), am = ap.a0 + ap.sweep / 2; mx = ap.cx + ap.r * Math.cos(am); my = ap.cy + ap.r * Math.sin(am); txt = "R" + +ap.r.toFixed(2); }
        else { mx = (e.ax + e.bx) / 2; my = (e.ay + e.by) / 2; txt = String(+len.toFixed(2)); }
        const q = w2s(mx, my), tw = ctx.measureText(txt).width + 6;
        ctx.fillStyle = "rgba(255,255,255,.85)"; ctx.fillRect(q[0] - tw / 2, q[1] - 8, tw, 16);
        ctx.fillStyle = "#33506e"; ctx.fillText(txt, q[0], q[1]);
      }
    }
    // étiquettes de matériau
    ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "left"; ctx.textBaseline = "middle";
    for (const l of S.labels) {
      const q = w2s(l.x, l.y), m = FT.getMaterial(S, l.mat);
      let txt = m.name.split(" (")[0];
      if (S.physics === "mag") {
        if (l.srcType === "NI") txt += "  " + l.srcVal + " A·t";
        if (l.srcType === "J") txt += "  " + l.srcVal + " A/mm²";
      } else {
        if (l.elecType === "V") txt += "  " + l.elecVal + " V";
        if (l.elecType === "rho") txt += "  ρ=" + l.elecVal;
      }
      ctx.beginPath(); ctx.arc(q[0], q[1], 4, 0, 6.3);
      ctx.fillStyle = result ? "#fff" : "#1463d6"; ctx.fill();
      ctx.strokeStyle = "#fff"; ctx.lineWidth = 3; ctx.strokeText(txt, q[0] + 8, q[1]);
      ctx.fillStyle = "#12202e"; ctx.fillText(txt, q[0] + 8, q[1]);
      if (m.Br > 0 && S.physics === "mag") {
        const an = ((+l.magAngle || 0) * Math.PI) / 180, dx = Math.cos(an) * 18, dy = -Math.sin(an) * 18;
        arrow(q[0] - dx, q[1] - dy + 16, q[0] + dx, q[1] + dy + 16, "#b91c1c", 2);
      }
    }
    // tracé en cours
    if (act && act.type === "draw" && act.pts.length > 1) {
      ctx.beginPath();
      act.pts.forEach((p, k) => { const q = w2s(p.x, p.y); k ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1]); });
      ctx.strokeStyle = "#1463d6"; ctx.lineWidth = 2.5; ctx.stroke();
    }
    if (act && act.type === "calib" && act.p1) {
      const q = w2s(act.p1.x, act.p1.y);
      ctx.beginPath(); ctx.arc(q[0], q[1], 6, 0, 6.3); ctx.strokeStyle = "#e11d48"; ctx.lineWidth = 2; ctx.stroke();
    }
    if (probe && result) {
      const q = w2s(probe.x, probe.y);
      ctx.beginPath(); ctx.moveTo(q[0] - 9, q[1]); ctx.lineTo(q[0] + 9, q[1]); ctx.moveTo(q[0], q[1] - 9); ctx.lineTo(q[0], q[1] + 9);
      ctx.strokeStyle = "#fff"; ctx.lineWidth = 3; ctx.stroke(); ctx.strokeStyle = "#000"; ctx.lineWidth = 1.2; ctx.stroke();
    }
    drawScaleBar();
    if (result && post) drawColorBar();
  }
  function arrow(x0, y0, x1, y1, col, w) {
    const an = Math.atan2(y1 - y0, x1 - x0), h = 3 + 2.2 * w;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1);
    ctx.moveTo(x1 - h * Math.cos(an - 0.5), y1 - h * Math.sin(an - 0.5)); ctx.lineTo(x1, y1); ctx.lineTo(x1 - h * Math.cos(an + 0.5), y1 - h * Math.sin(an + 0.5));
    ctx.strokeStyle = col; ctx.lineWidth = w; ctx.stroke();
  }
  function fieldAt(x, y) {
    const bary = [0, 0, 0], t = mesh.locator.find(x, y, bary);
    if (t < 0) return null;
    const ip = (s) => bary[0] * s[3 * t] + bary[1] * s[3 * t + 1] + bary[2] * s[3 * t + 2];
    return { t, vx: ip(post.sx), vy: ip(post.sy), bary };
  }
  function drawArrows() {
    const step = 34, ref = result.type === "mag" ? post.refV : post.refV;
    for (let sy = step / 2; sy < H; sy += step)
      for (let sx = step / 2; sx < W; sx += step) {
        const p = s2w(sx, sy), f = fieldAt(p.x, p.y);
        if (!f) continue;
        const m = Math.hypot(f.vx, f.vy);
        if (!(m > 0.03 * ref)) continue;
        const L = 7 + 7 * Math.min(1, m / ref), ux = (f.vx / m) * L, uy = (-f.vy / m) * L;
        arrow(sx - ux, sy - uy, sx + ux, sy + uy, "rgba(0,0,0,.8)", 1.2);
      }
  }
  function drawScaleBar() {
    const len = niceStep(110 / view.z), px = len * view.z, x = 14, y = H - 16;
    ctx.beginPath(); ctx.moveTo(x, y - 5); ctx.lineTo(x, y); ctx.lineTo(x + px, y); ctx.lineTo(x + px, y - 5);
    ctx.strokeStyle = "#fff"; ctx.lineWidth = 4; ctx.stroke(); ctx.strokeStyle = "#1c2733"; ctx.lineWidth = 1.5; ctx.stroke();
    ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "center"; ctx.textBaseline = "bottom";
    ctx.strokeStyle = "#fff"; ctx.lineWidth = 3; ctx.strokeText(mm(len), x + px / 2, y - 4);
    ctx.fillStyle = "#1c2733"; ctx.fillText(mm(len), x + px / 2, y - 4);
  }
  function drawColorBar() {
    const h = Math.min(220, H - 90), x = W - 78, y = 34, c = [0, 0, 0, 0];
    ctx.fillStyle = "rgba(255,255,255,.92)"; ctx.fillRect(x - 10, y - 26, 82, h + 40);
    for (let i = 0; i < h; i++) { cmap(1 - i / (h - 1), c, 0); ctx.fillStyle = "rgb(" + (c[0] | 0) + "," + (c[1] | 0) + "," + (c[2] | 0) + ")"; ctx.fillRect(x, y + i, 14, 1.5); }
    ctx.fillStyle = "#1c2733"; ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "left"; ctx.textBaseline = "middle";
    for (let k = 0; k <= 4; k++) { const v = post.hi - ((post.hi - post.lo) * k) / 4; ctx.fillText(+v.toPrecision(3) + "", x + 19, y + (h * k) / 4); }
    ctx.font = "600 12px system-ui, sans-serif"; ctx.fillText(post.title + " (" + post.unit + ")", x - 4, y - 13);
  }

  // ---------- modèle : modifications, annulation, sauvegarde ----------
  const snapshot = () => JSON.stringify({ edges: S.edges, labels: S.labels });
  function pushUndo() { undo.push(snapshot()); if (undo.length > 80) undo.shift(); redo.length = 0; }
  function restore(js) { const o = JSON.parse(js); S.edges = o.edges; S.labels = o.labels; sel = -1; changed(); }
  function save() {
    try { localStorage.setItem(KEY, JSON.stringify({ S, view })); } catch (e) { /* stockage indisponible */ }
  }
  function changed(fast) {
    result = null; post = null; probe = null; $("readout").hidden = true;
    save();
    clearTimeout(meshTimer);
    meshTimer = setTimeout(remesh, fast ? 0 : 160);
    updateUI();
    draw();
  }
  function remesh() {
    clearTimeout(meshTimer); meshTimer = 0;
    meshError = "";
    try { mesh = FT.buildMesh(S); } catch (e) { mesh = null; meshError = e.message; console.error(e); }
    uploadMesh(); uploadColors();
    updateUI();
    draw();
  }
  function updateUI() {
    $("btnDel").disabled = sel < 0;
    $("btnUndo").disabled = !undo.length; $("btnRedo").disabled = !redo.length;
    $("btnCalib").disabled = $("btnImgDel").disabled = !bg;
    $("meshInfo").textContent = meshError ? "Erreur de maillage : " + meshError
      : mesh ? mesh.x.length + " nœuds, " + mesh.tri.length / 3 + " éléments, " + mesh.nReg + " zones" + (mesh.warnings.length ? " — " + mesh.warnings.join(" ; ") : "")
      : "Aucun maillage";
    if (!result) $("results").textContent = "Dessinez un circuit, affectez les matériaux (outil Matériau) puis lancez le calcul.";
    let msg = "";
    if (sel >= 0 && S.edges[sel]) {
      const e = S.edges[sel];
      msg = e.type === "arc" ? "Arc : rayon " + mm(FT.arcParams(e).r) + ", longueur " + mm(FT.edgeLength(e)) : "Segment : longueur " + mm(FT.edgeLength(e));
    }
    $("stMsg").textContent = msg || { draw: S.drawMode === "free" ? "Mode Patate : le tracé est gardé tel quel ; ses extrémités s'accrochent aux traits voisins." : "Mode Grille : droites, arcs, cercles et contours s'accrochent à la grille.", edit: "Touchez une arête pour la sélectionner, faites glisser un sommet pour le déplacer.", erase: "Passez sur les arêtes à effacer.", mat: "Touchez une zone fermée pour choisir son matériau et ses sources.", probe: "Touchez le résultat pour lire les valeurs locales.", pan: "Faites glisser pour déplacer, molette ou pincement pour zoomer.", calib: "Touchez deux points de l'image dont vous connaissez la distance." }[tool];
  }
  function setTool(t) {
    tool = t; act = null;
    document.querySelectorAll("#tools button").forEach((b) => b.classList.toggle("on", b.dataset.tool === t));
    cv.style.cursor = t === "pan" ? "grab" : "crosshair";
    if (t !== "edit") sel = -1;
    // retour à l'édition : on retire l'affichage du résultat
    if (result && (t === "draw" || t === "edit" || t === "erase")) {
      result = null; post = null; probe = null; $("readout").hidden = true;
      uploadColors();
    }
    updateUI(); draw();
  }
  function fit() {
    if (!S.edges.length) { view = { x: 30, y: 30, z: 6 }; draw(); return; }
    const bb = FT.modelBBox(S.edges), w = bb.x1 - bb.x0 || 1, h = bb.y1 - bb.y0 || 1;
    view.x = (bb.x0 + bb.x1) / 2; view.y = (bb.y0 + bb.y1) / 2;
    view.z = Math.min(W / (w * 1.5), H / (h * 1.5));
    draw();
  }

  // ---------- calcul ----------
  function solve() {
    if (!S.edges.length) { toast("Dessinez d'abord une géométrie."); return; }
    toast("Calcul en cours…", 60000);
    setTimeout(() => {
      try {
        if (meshTimer || !mesh) remesh();
        if (!mesh) { toast("Maillage impossible : " + (meshError || "géométrie vide")); return; }
        result = FT.solve(S, mesh);
        result.props = FT.regionProps(S, mesh, S.physics);
        buildPost();
        let ref = 0;
        const vx = post.sx, vy = post.sy;
        for (let k = 0; k < vx.length; k += 7) ref = Math.max(ref, Math.hypot(vx[k], vy[k]));
        post.refV = ref;
        uploadColors();
        showResults();
        $("toast").hidden = true;
      } catch (e) {
        console.error(e); result = null; post = null;
        toast("Erreur de calcul : " + e.message, 6000);
      }
      draw();
    }, 30);
  }
  function showResults() {
    const r = result, d = S.depth * 1e-3, out = [];
    out.push("<b>Maillage</b> : " + mesh.x.length + " nœuds, " + mesh.tri.length / 3 + " éléments");
    if (r.type === "mag") {
      out.push("<b>Solveur</b> : " + (r.nonlinear ? "non linéaire, " + r.iters + " itérations de Newton" : "linéaire") + ", " + r.ms + " ms");
      if (!r.converged) out.push('<span class="warn">Convergence incomplète (résidu ' + r.relRes.toExponential(1) + ")</span>");
      let pos = 0, neg = 0, hasPM = false;
      result.props.forEach((p, i) => { const I = p.J * mesh.regArea[i]; if (I > 0) pos += I; else neg -= I; if (p.brx || p.bry) hasPM = true; });
      const NI = Math.max(pos, neg);
      if (!NI && !hasPM) out.push('<span class="warn">Aucune source : ajoutez des ampères-tours ou un aimant avec l\'outil Matériau.</span>');
      if (NI && Math.abs(pos - neg) > 1e-6 * NI) out.push('<span class="warn">Somme des courants non nulle (' + fmt(pos - neg, "A·t") + ") : le retour se fait par la frontière.</span>");
      out.push("<b>|B| max affiché</b> : " + fmt(post.hi, post.unit));
      out.push("<b>Énergie magnétique</b> : " + fmt(r.energy * d, "J") + " (profondeur " + S.depth + " mm)");
      if (NI) {
        const phi = (r.linkage * d) / NI;
        out.push("<b>Ampères-tours</b> : " + fmt(NI, "A·t"));
        out.push("<b>Flux lié par spire</b> : " + fmt(phi, "Wb"));
        out.push("<b>A<sub>L</sub> = Φ/NI</b> : " + fmt(phi / NI, "H/spire²") + (r.nonlinear || hasPM ? " (apparente)" : ""));
      }
    } else {
      out.push("<b>Solveur</b> : " + r.cgIt + " itérations, " + r.ms + " ms");
      if (r.pinned) out.push('<span class="warn">Aucun potentiel imposé : ajoutez-en avec l\'outil Matériau.</span>');
      let vmin = Infinity, vmax = -Infinity;
      result.props.forEach((p) => { if (p.V != null) { vmin = Math.min(vmin, p.V); vmax = Math.max(vmax, p.V); } });
      if (S.bcElec === "V0") { vmin = Math.min(vmin, 0); vmax = Math.max(vmax, 0); }
      out.push("<b>Énergie électrique</b> : " + fmt(r.energy * d, "J") + " (profondeur " + S.depth + " mm)");
      if (vmax > vmin) out.push("<b>Capacité 2W/ΔV²</b> : " + fmt((2 * r.energy * d) / ((vmax - vmin) * (vmax - vmin)), "F") + " (ΔV = " + (vmax - vmin) + " V)");
    }
    out.push("");
    out.push("Outil <b>Sonde</b> : touchez une zone pour ses valeurs locales, sa moyenne et la force qu'elle subit.");
    $("results").innerHTML = out.join("\n");
  }
  function updateProbe(p) {
    if (!result || !post) { toast("Lancez d'abord le calcul."); return; }
    const f = fieldAt(p.x, p.y);
    if (!f) return;
    probe = p;
    const r = result, t = f.t, reg = mesh.reg[t], lab = mesh.regLabel[reg];
    const L = ["x = " + p.x.toFixed(2) + " mm   y = " + p.y.toFixed(2) + " mm"];
    if (r.type === "mag") {
      const A = f.bary[0] * r.A[mesh.tri[3 * t]] + f.bary[1] * r.A[mesh.tri[3 * t + 1]] + f.bary[2] * r.A[mesh.tri[3 * t + 2]];
      L.push("|B| = " + fmt(Math.hypot(f.vx, f.vy), "T"), "Bx = " + fmt(f.vx, "T") + "   By = " + fmt(f.vy, "T"),
        "|H| = " + fmt(r.Hm[t], "A/m") + "   µr = " + +(1 / (r.nu[t] * FT.MU0)).toPrecision(4), "A = " + fmt(A, "Wb/m"));
    } else {
      const V = f.bary[0] * r.V[mesh.tri[3 * t]] + f.bary[1] * r.V[mesh.tri[3 * t + 1]] + f.bary[2] * r.V[mesh.tri[3 * t + 2]];
      L.push("V = " + fmt(V, "V"), "|E| = " + fmt(Math.hypot(f.vx, f.vy) * 1e-3, "V/mm"), "Ex = " + fmt(f.vx * 1e-3, "V/mm") + "   Ey = " + fmt(f.vy * 1e-3, "V/mm"), "|D| = " + fmt(r.Dm[t], "C/m²"));
    }
    const st = post.stats[reg] || (post.stats[reg] = FT.regionStats(mesh, r, reg, S.depth * 1e-3));
    L.push("", "Zone : " + FT.getMaterial(S, lab ? lab.mat : "air").name, "Aire = " + +st.area.toPrecision(4) + " mm²",
      (r.type === "mag" ? "|B| moyen = " + fmt(st.mean, "T") : "|E| moyen = " + fmt(st.mean * 1e-3, "V/mm")),
      "Énergie = " + fmt(st.energy, "J"));
    if (lab) L.push("Force : Fx = " + fmt(st.fx, "N") + "   Fy = " + fmt(st.fy, "N") + (st.contact ? "\n(zone non entourée d'air : force indicative)" : ""));
    const ro = $("readout");
    ro.textContent = L.join("\n"); ro.hidden = false;
    draw();
  }

  // ---------- boîte de dialogue matériau ----------
  let editing = null;
  function fillMatSelect(selId) {
    const s = $("dMat");
    s.innerHTML = "";
    for (const m of FT.defaultMaterials.concat(S.customMats)) { const o = document.createElement("option"); o.value = m.id; o.textContent = m.name; s.appendChild(o); }
    s.value = selId;
  }
  function openRegion(p) {
    if (meshTimer) remesh();
    if (!mesh) { toast("Dessinez d'abord un contour fermé."); return; }
    const t = mesh.locator.find(p.x, p.y);
    if (t < 0) { toast("Point hors du domaine."); return; }
    const lab = mesh.regLabel[mesh.reg[t]];
    editing = { lab, p };
    const l = lab || { mat: "steel", srcType: "none", srcVal: 0, magAngle: 0, elecType: "none", elecVal: 0 };
    fillMatSelect(l.mat);
    $("dAngle").value = l.magAngle || 0; $("dSrc").value = l.srcType || "none"; $("dSrcVal").value = l.srcVal || 0;
    $("dElec").value = l.elecType || "none"; $("dElecVal").value = l.elecVal || 0;
    $("dReset").hidden = !lab;
    $("dCustom").open = false;
    syncDlg();
    $("dlg").showModal();
  }
  // pose du matériau de la palette ; toucher à nouveau une zone déjà de ce matériau ouvre ses réglages
  function paintRegion(p) {
    if (meshTimer) remesh();
    if (!mesh) { toast("Dessinez d'abord un contour fermé."); return; }
    const t = mesh.locator.find(p.x, p.y);
    if (t < 0) { toast("Point hors du domaine."); return; }
    const lab = mesh.regLabel[mesh.reg[t]], m = FT.getMaterial(S, armed);
    if (lab && lab.mat === armed) { openRegion(p); return; }
    if (armed === "air" && !lab) return;
    pushUndo();
    if (armed === "air") S.labels = S.labels.filter((l) => l !== lab);
    else if (lab) lab.mat = armed;
    else S.labels.push({ x: p.x, y: p.y, mat: armed, srcType: "none", srcVal: 0, magAngle: 0, elecType: "none", elecVal: 0 });
    changed(true);
    toast(m.name + (armed === "air" ? "" : " — touchez à nouveau pour les sources"), 1800);
  }
  const shortName = (n) => n.replace(/\s*\(.*\)\s*$/, "");
  function renderPalette() {
    const box = $("palList");
    box.innerHTML = "";
    const add = (m) => {
      const b = document.createElement("button");
      b.type = "button"; b.className = armed === m.id ? "on" : "";
      b.title = m.name + " — µr " + m.mur + (m.Js ? ", Js " + m.Js + " T" : "") + (m.Br ? ", Br " + m.Br + " T" : "") + ", εr " + m.epsr + (m.cond ? ", conducteur" : "");
      b.innerHTML = '<span class="dot"></span><span></span>';
      b.firstChild.style.background = m.color; b.lastChild.textContent = shortName(m.name);
      b.onclick = () => {
        armed = armed === m.id ? null : m.id;
        renderPalette();
        if (armed) { if (tool !== "mat") setTool("mat"); toast("Touchez une zone pour y mettre : " + m.name, 1800); }
      };
      box.appendChild(b);
    };
    FT.defaultMaterials.forEach(add);
    if (S.customMats.length) { const sp = document.createElement("div"); sp.className = "sep"; box.appendChild(sp); S.customMats.forEach(add); }
  }
  function setPalette(show) {
    $("palette").classList.toggle("hidden", !show); $("palShow").hidden = show;
    try { localStorage.setItem("fieldtool.palette", show ? "1" : "0"); } catch (e) { /* stockage indisponible */ }
    resize();
  }
  $("palHide").onclick = () => setPalette(false);
  $("palShow").onclick = () => setPalette(true);
  $("palAdd").onclick = () => {
    $("nName").value = "";
    $("nColor").value = ["#c7b8ea", "#a8d5ba", "#f2c6a0", "#9fc5e8", "#e6b8af", "#d5a6bd"][S.customMats.length % 6];
    $("dNew").showModal();
  };
  $("nCancel").onclick = () => $("dNew").close();
  $("nOk").onclick = () => {
    const m = { id: "c" + Date.now().toString(36), name: $("nName").value.trim() || "Matériau " + (S.customMats.length + 1), color: $("nColor").value,
      mur: Math.max(1e-3, num($("nMur").value, 1)), Js: Math.max(0, num($("nJs").value, 0)), Br: Math.max(0, num($("nBr").value, 0)),
      epsr: Math.max(1e-3, num($("nEps").value, 1)), cond: $("nCond").checked };
    if (m.Br > 0) m.Js = 0;
    S.customMats.push(m); save();
    $("dNew").close();
    armed = m.id; renderPalette();
    if (tool !== "mat") setTool("mat");
    toast("« " + m.name + " » créé : touchez une zone pour l'y mettre.", 2500);
  };
  function syncDlg() {
    $("dMagnet").hidden = !(FT.getMaterial(S, $("dMat").value).Br > 0) || S.physics !== "mag";
    $("dFsMag").hidden = S.physics !== "mag";
    $("dFsElec").hidden = S.physics !== "elec";
  }
  function closeDlg(apply) {
    if (apply && editing) {
      pushUndo();
      const v = { mat: $("dMat").value, magAngle: num($("dAngle").value, 0), srcType: $("dSrc").value, srcVal: num($("dSrcVal").value, 0), elecType: $("dElec").value, elecVal: num($("dElecVal").value, 0) };
      if (apply === "reset") S.labels = S.labels.filter((l) => l !== editing.lab);
      else if (editing.lab) Object.assign(editing.lab, v);
      else S.labels.push(Object.assign({ x: editing.p.x, y: editing.p.y }, v));
      changed(true);
    }
    editing = null;
    $("dlg").close();
  }

  // ---------- entrées pointeur ----------
  const touches = new Map();
  let pinch = null;
  const evPt = (e) => { const r = cv.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  function hitVertex(pt) {
    let best = null, bd = 16;
    for (const e of S.edges)
      for (const v of [[e.ax, e.ay], [e.bx, e.by]]) {
        const q = w2s(v[0], v[1]), d = Math.hypot(q[0] - pt.x, q[1] - pt.y);
        if (d < bd) { bd = d; best = { x: v[0], y: v[1] }; }
      }
    return best;
  }
  function hitEdge(pt, tolPx) {
    const p = s2w(pt.x, pt.y);
    let best = -1, bd = tolPx / view.z;
    S.edges.forEach((e, i) => { const d = FT.distToEdge(e, p.x, p.y); if (d < bd) { bd = d; best = i; } });
    return best;
  }
  cv.addEventListener("pointerdown", (e) => {
    e.preventDefault();
    try { cv.setPointerCapture(e.pointerId); } catch (err) { /* pointeur déjà relâché */ }
    const pt = evPt(e);
    if (e.pointerType === "pen") {
      penSeen = true;
      // bouton latéral (2) ou gomme (5) d'un stylet : bascule Grille / Patate
      if (e.button === 2 || e.button === 5) { setDrawMode(S.drawMode === "free" ? "grid" : "free"); return; }
    }
    if (e.pointerType === "touch") {
      touches.set(e.pointerId, pt);
      if (touches.size === 2) {
        act = null;
        const [a, b] = [...touches.values()];
        pinch = { d: Math.hypot(a.x - b.x, a.y - b.y), cx: (a.x + b.x) / 2, cy: (a.y + b.y) / 2 };
        // mémorisé pour reconnaître un tap à deux doigts (bascule Grille / Patate)
        pinch.tap = { t0: performance.now(), d0: pinch.d, cx0: pinch.cx, cy0: pinch.cy, z0: view.z, x0: view.x, y0: view.y, ok: true };
        draw();
        return;
      }
      if (touches.size > 2) return;
      if (penSeen && tool !== "mat" && tool !== "probe") { act = { type: "pan", last: pt }; return; }
    }
    if ((e.pointerType === "mouse" && e.button !== 0) || tool === "pan") { act = { type: "pan", last: pt }; return; }
    const w = s2w(pt.x, pt.y);
    if (tool === "draw") act = { type: "draw", pts: [w] };
    else if (tool === "erase") { act = { type: "erase", hit: new Set() }; eraseAt(pt); }
    else if (tool === "mat") act = { type: "tap", start: pt, pts: [w], touch: e.pointerType === "touch" };
    else if (tool === "probe") { act = { type: "probe" }; updateProbe(w); }
    else if (tool === "calib") act = { type: "calib", p1: act && act.p1, start: pt };
    else if (tool === "edit") {
      const v = hitVertex(pt);
      if (v) act = { type: "drag", from: v, cur: v, before: snapshot() };
      else { sel = hitEdge(pt, 14); updateUI(); }
    }
    draw();
  });
  function eraseAt(pt) { const i = hitEdge(pt, 12); if (i >= 0) act.hit.add(i); }
  cv.addEventListener("pointermove", (e) => {
    const pt = evPt(e);
    const w = s2w(pt.x, pt.y), sp = snapPt(w, act && act.type === "drag" ? act.cur : null);
    $("stPos").textContent = "x " + sp.x.toFixed(2) + "  y " + sp.y.toFixed(2) + " mm";
    if (e.pointerType === "touch" && touches.has(e.pointerId)) {
      touches.set(e.pointerId, pt);
      if (pinch && touches.size === 2) {
        const [a, b] = [...touches.values()];
        const d = Math.hypot(a.x - b.x, a.y - b.y), cx = (a.x + b.x) / 2, cy = (a.y + b.y) / 2;
        const anchor = s2w(pinch.cx, pinch.cy);
        view.z = Math.min(4000, Math.max(0.02, (view.z * d) / (pinch.d || d)));
        view.x = anchor.x - (cx - W / 2) / view.z; view.y = anchor.y + (cy - H / 2) / view.z;
        const tap = pinch.tap;
        if (tap && (Math.abs(d - tap.d0) > 14 || Math.hypot(cx - tap.cx0, cy - tap.cy0) > 14)) tap.ok = false;
        pinch = { d, cx, cy, tap };
        draw();
        return;
      }
    }
    if (!act) return;
    if (act.type === "pan") { view.x -= (pt.x - act.last.x) / view.z; view.y += (pt.y - act.last.y) / view.z; act.last = pt; }
    else if (act.type === "draw") {
      const ev = e.getCoalescedEvents ? e.getCoalescedEvents() : [];
      if (ev.length > 1) for (const c of ev) { const q = evPt(c); act.pts.push(s2w(q.x, q.y)); }
      else act.pts.push(w);
    } else if (act.type === "tap") {
      act.pts.push(w);
      if (Math.hypot(pt.x - act.start.x, pt.y - act.start.y) > 12) {
        // un doigt après le stylet déplace la vue ; un trait de stylet ou de souris redevient un tracé
        if (act.touch && penSeen) act = { type: "pan", last: pt };
        else { act.type = "draw"; act.fromMat = true; }
      }
    } else if (act.type === "erase") eraseAt(pt);
    else if (act.type === "probe") updateProbe(w);
    else if (act.type === "drag") {
      if (sp.x !== act.cur.x || sp.y !== act.cur.y) {
        for (const ed of S.edges) {
          if (ed.ax === act.cur.x && ed.ay === act.cur.y) { ed.ax = sp.x; ed.ay = sp.y; }
          if (ed.bx === act.cur.x && ed.by === act.cur.y) { ed.bx = sp.x; ed.by = sp.y; }
        }
        act.cur = sp; act.moved = true;
      }
    }
    draw();
  });
  function endPointer(e) {
    const pt = evPt(e);
    if (e.pointerType === "touch") { touches.delete(e.pointerId); if (touches.size < 2 && pinch) {
        const tap = pinch.tap;
        pinch = null;
        if (tap && tap.ok && performance.now() - tap.t0 < 350) {
          // tap à deux doigts : on annule le petit mouvement de vue éventuel et on bascule le mode
          view.z = tap.z0; view.x = tap.x0; view.y = tap.y0;
          setDrawMode(S.drawMode === "free" ? "grid" : "free");
        }
        save(); return;
      } }
    const a = act;
    if (!a) return;
    if (a.type === "calib") { act = { type: "calib", p1: a.p1 }; calibTap(s2w(pt.x, pt.y)); return; }
    act = null;
    if (a.type === "draw") {
      const ne = S.drawMode === "free"
        ? FT.freehandStroke(a.pts, { tol: 1.2 / view.z, snapR: 12 / view.z, edges: S.edges })
        : FT.recognizeStroke(a.pts, S.grid);
      if (ne.length) { pushUndo(); S.edges.push(...ne); changed(); }
      if (a.fromMat) { setTool("draw"); toast("Retour au dessin", 1200); }
      else if (a.pts.length > 6) toast(S.drawMode === "free" ? "Tracé trop petit." : "Tracé trop petit par rapport au pas de grille (" + S.grid + " mm).");
    } else if (a.type === "erase") {
      if (a.hit.size) { pushUndo(); S.edges = S.edges.filter((_, i) => !a.hit.has(i)); changed(); }
    } else if (a.type === "tap") {
      if (Math.hypot(pt.x - a.start.x, pt.y - a.start.y) < 12) (armed ? paintRegion : openRegion)(s2w(pt.x, pt.y));
    } else if (a.type === "drag" && a.moved) {
      undo.push(a.before); redo.length = 0;
      S.edges = S.edges.filter((ed) => ed.ax !== ed.bx || ed.ay !== ed.by);
      changed();
    } else if (a.type === "pan") save();
    draw();
  }
  cv.addEventListener("pointerup", endPointer);
  cv.addEventListener("pointercancel", (e) => { touches.delete(e.pointerId); pinch = null; if (act && act.type !== "calib") act = null; draw(); });
  cv.addEventListener("contextmenu", (e) => e.preventDefault());
  cv.addEventListener("wheel", (e) => {
    e.preventDefault();
    const pt = evPt(e), anchor = s2w(pt.x, pt.y);
    view.z = Math.min(4000, Math.max(0.02, view.z * Math.exp(-e.deltaY * 0.0015)));
    view.x = anchor.x - (pt.x - W / 2) / view.z; view.y = anchor.y + (pt.y - H / 2) / view.z;
    draw();
  }, { passive: false });

  // ---------- image de fond ----------
  function calibTap(p) {
    if (!act.p1) { act.p1 = p; toast("Touchez le second point."); draw(); return; }
    const d0 = Math.hypot(p.x - act.p1.x, p.y - act.p1.y), p1 = act.p1;
    setTool("draw");
    if (d0 <= 0) return;
    ask({ msg: "Distance réelle entre les deux points (mm) :", input: 10, buttons: [["Annuler", null], ["Valider", true, 1]] }).then((v) => {
      const real = num(v, 0);
      if (real > 0 && bg) {
        const k = real / d0;
        bg.mm *= k; bg.x = p1.x + (bg.x - p1.x) * k; bg.y = p1.y + (bg.y - p1.y) * k;
        toast("Image mise à l'échelle : " + (1 / bg.mm).toFixed(2) + " px/mm.");
      }
      draw();
    });
    draw();
  }
  $("btnImg").onclick = () => $("fileImg").click();
  $("fileImg").onchange = (e) => {
    const f = e.target.files[0];
    if (!f) return;
    const img = new Image();
    img.onload = () => {
      const mmpp = (0.7 * W) / view.z / img.width;
      bg = { img, mm: mmpp, x: view.x - (img.width * mmpp) / 2, y: view.y + (img.height * mmpp) / 2, alpha: num($("inAlpha").value, 0.5) };
      updateUI(); draw();
      toast("Image importée. Utilisez « Mettre à l'échelle » puis tracez par-dessus.", 4000);
    };
    img.src = URL.createObjectURL(f);
    e.target.value = "";
  };
  $("btnCalib").onclick = () => { setTool("calib"); act = { type: "calib", p1: null }; };
  $("btnImgDel").onclick = () => { bg = null; updateUI(); draw(); };
  $("inAlpha").oninput = (e) => { if (bg) { bg.alpha = num(e.target.value, 0.5); draw(); } };

  // ---------- commandes ----------
  document.querySelectorAll("#tools button").forEach((b) => (b.onclick = () => setTool(b.dataset.tool)));
  function setDrawMode(m, quiet) {
    S.drawMode = m;
    document.querySelectorAll("#dmode button").forEach((b) => b.classList.toggle("on", b.dataset.dm === m));
    if (!quiet) { save(); toast(m === "free" ? "Mode Patate : tracé libre" : "Mode Grille : tracé accroché à la grille", 1500); }
    updateUI(); draw();
  }
  document.querySelectorAll("#dmode button").forEach((b) => (b.onclick = () => { setDrawMode(b.dataset.dm); if (tool !== "draw") setTool("draw"); }));
  function setPhysics(p) {
    S.physics = p;
    document.querySelectorAll("#phys button").forEach((b) => b.classList.toggle("on", b.dataset.phys === p));
    const q = $("inQty");
    q.innerHTML = "";
    for (const [v, t] of QTY[p]) { const o = document.createElement("option"); o.value = v; o.textContent = t; q.appendChild(o); }
    disp.qty = QTY[p][0][0]; disp.max = null; $("inMax").value = "";
    $("ckLines").parentNode.lastChild.textContent = p === "mag" ? " Lignes de flux" : " Équipotentielles";
  }
  document.querySelectorAll("#phys button").forEach((b) => (b.onclick = () => { setPhysics(b.dataset.phys); result = null; post = null; uploadColors(); changed(true); }));
  $("btnSolve").onclick = solve;
  $("btnUndo").onclick = () => { if (undo.length) { redo.push(snapshot()); restore(undo.pop()); } };
  $("btnRedo").onclick = () => { if (redo.length) { undo.push(snapshot()); restore(redo.pop()); } };
  $("btnDel").onclick = () => { if (sel >= 0) { pushUndo(); S.edges.splice(sel, 1); sel = -1; changed(); } };
  $("btnFit").onclick = fit;
  $("btnPanel").onclick = () => { $("panel").classList.toggle("hidden"); resize(); };
  $("inGrid").onchange = (e) => { S.grid = Math.max(1e-3, num(e.target.value, S.grid)); e.target.value = S.grid; save(); draw(); };
  $("inGridReal").onchange = (e) => { S.gridReal = Math.max(1e-3, num(e.target.value, S.gridReal)); e.target.value = S.gridReal; save(); draw(); };
  $("inDepth").onchange = (e) => { S.depth = Math.max(1e-3, num(e.target.value, S.depth)); e.target.value = S.depth; save(); if (result) { post.stats = {}; showResults(); } };
  $("inMargin").onchange = (e) => { S.airMargin = Math.min(10, Math.max(0.2, num(e.target.value, 1))); e.target.value = S.airMargin; changed(); };
  $("inMesh").onchange = (e) => { S.meshN = num(e.target.value, 40); changed(); };
  $("inBc").onchange = (e) => { S.bcElec = e.target.value; changed(true); };
  $("inQty").onchange = (e) => { disp.qty = e.target.value; disp.max = null; $("inMax").value = ""; refreshPost(); };
  $("inMax").onchange = (e) => { disp.max = num(e.target.value, null); refreshPost(); };
  function refreshPost() {
    if (!result) return;
    const ref = post && post.refV, st = post && post.stats;
    buildPost(); post.refV = ref; post.stats = st || {};
    uploadColors(); draw();
  }
  for (const [id, k] of [["ckLines", "lines"], ["ckArrows", "arrows"], ["ckMesh", "mesh"], ["ckDims", "dims"]]) $(id).onchange = (e) => { disp[k] = e.target.checked; draw(); };

  $("dMat").onchange = syncDlg;
  $("dOk").onclick = () => closeDlg(true);
  $("dCancel").onclick = () => closeDlg(false);
  $("dReset").onclick = () => closeDlg("reset");
  document.querySelectorAll("#dMagnet [data-ang]").forEach((b) => (b.onclick = () => ($("dAngle").value = b.dataset.ang)));
  $("cAdd").onclick = () => {
    const pal = ["#c7b8ea", "#a8d5ba", "#f2c6a0", "#9fc5e8", "#e6b8af", "#d5a6bd"];
    const m = { id: "c" + Date.now().toString(36), name: $("cName").value.trim() || "Matériau " + (S.customMats.length + 1), color: pal[S.customMats.length % pal.length],
      mur: Math.max(1e-3, num($("cMur").value, 1)), Js: Math.max(0, num($("cJs").value, 0)), Br: Math.max(0, num($("cBr").value, 0)), epsr: Math.max(1e-3, num($("cEps").value, 1)), cond: $("cCond").checked };
    if (m.Br > 0) m.Js = 0;
    S.customMats.push(m);
    fillMatSelect(m.id); syncDlg(); renderPalette();
    $("dCustom").open = false;
  };

  function loadProject(o) {
    S = Object.assign(defaults(), o.S || o);
    if (o.view) view = o.view;
    undo.length = redo.length = 0; sel = -1; mesh = null;
    syncInputs();
    changed(true);
  }
  function syncInputs() {
    $("inGrid").value = S.grid; $("inGridReal").value = S.gridReal; $("inDepth").value = S.depth; $("inMargin").value = S.airMargin; $("inMesh").value = S.meshN; $("inBc").value = S.bcElec;
    setPhysics(S.physics);
    setDrawMode(S.drawMode || "grid", true);
    if (armed && !FT.defaultMaterials.concat(S.customMats).some((m) => m.id === armed)) armed = null;
    renderPalette();
  }
  // ---------- exemples ----------
  function exampleKit() {
    const edges = [];
    const poly = (p) => p.forEach((a, i) => { const b = p[(i + 1) % p.length]; edges.push({ type: "line", ax: a[0], ay: a[1], bx: b[0], by: b[1] }); });
    const rect = (x0, y0, x1, y1) => poly([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]);
    const circle = (cx, cy, r) => {
      edges.push({ type: "arc", ax: cx + r, ay: cy, bx: cx - r, by: cy, bulge: 1 });
      edges.push({ type: "arc", ax: cx - r, ay: cy, bx: cx + r, by: cy, bulge: 1 });
    };
    const L = (x, y, mat, o = {}) => ({
      x, y, mat, srcType: o.NI ? "NI" : o.J ? "J" : "none", srcVal: o.NI || o.J || 0, magAngle: o.ang || 0,
      elecType: o.V != null ? "V" : "none", elecVal: o.V || 0,
    });
    const make = (labels, opts) => ({ S: Object.assign(defaults(), opts, { edges, labels }) });
    return { poly, rect, circle, L, make };
  }
  const EXAMPLES = [
    {
      name: "Électroaimant en C", desc: "Acier doux, entrefer 4 mm, bobine 800 A·t",
      build() {
        const { poly, L, make } = exampleKit();
        poly([[0, 0], [60, 0], [60, 28], [45, 28], [45, 15], [15, 15], [15, 45], [45, 45], [45, 32], [60, 32], [60, 60], [0, 60]]);
        poly([[15, 20], [22, 20], [22, 40], [15, 40]]);
        poly([[-7, 20], [0, 20], [0, 40], [-7, 40]]);
        return make([L(7, 7, "steel"), L(18.5, 30, "cu", { NI: 800 }), L(-3.5, 30, "cu", { NI: -800 })], { grid: 5 });
      },
    },
    {
      name: "Inductance ferrite E-I", desc: "Entrefer 0,5 mm sur la jambe centrale, 50 A·t : lire A_L",
      build() {
        const { poly, rect, L, make } = exampleKit();
        poly([[0, 6], [5, 6], [5, 25], [15, 25], [15, 6.5], [25, 6.5], [25, 25], [35, 25], [35, 6], [40, 6], [40, 30], [0, 30]]);
        rect(0, 0, 40, 6);
        rect(7, 8.5, 14, 23);
        rect(26, 8.5, 33, 23);
        return make([L(2.5, 28, "ferrite"), L(20, 3, "ferrite"), L(10.5, 15, "cu", { NI: 50 }), L(29.5, 15, "cu", { NI: -50 })], { grid: 0.5, gridReal: 5, depth: 15 });
      },
    },
    {
      name: "Électroaimant de levage", desc: "Culasse en U et armature à 2 mm : force d'attraction",
      build() {
        const { poly, rect, L, make } = exampleKit();
        poly([[0, 10], [15, 10], [15, 35], [45, 35], [45, 10], [60, 10], [60, 50], [0, 50]]);
        rect(0, 0, 60, 8);
        rect(15, 14, 22, 33); rect(38, 14, 45, 33);
        rect(-7, 14, 0, 33); rect(60, 14, 67, 33);
        return make([
          L(30, 45, "steel"), L(30, 4, "steel"),
          L(18.5, 24, "cu", { NI: 1000 }), L(41.5, 24, "cu", { NI: 1000 }),
          L(-3.5, 24, "cu", { NI: -1000 }), L(63.5, 24, "cu", { NI: -1000 }),
        ], { grid: 1, gridReal: 10, depth: 30 });
      },
    },
    {
      name: "Circuit à aimant NdFeB", desc: "Aimant dans la culasse, entrefer 4 mm, sans bobine",
      build() {
        const { poly, rect, L, make } = exampleKit();
        poly([[0, 0], [50, 0], [50, 18], [40, 18], [40, 10], [10, 10], [10, 30], [40, 30], [40, 22], [50, 22], [50, 40], [0, 40]]);
        rect(0, 15, 10, 25);
        return make([L(25, 5, "steel"), L(25, 35, "steel"), L(5, 20, "ndfeb", { ang: 90 })], { grid: 1, gridReal: 10, depth: 20 });
      },
    },
    {
      name: "Condensateur plan FR4", desc: "Armatures 100 V / 0 V, 4 mm de FR4 : capacité et effets de bord",
      build() {
        const { rect, L, make } = exampleKit();
        rect(0, 6, 40, 8); rect(0, 2, 40, 6); rect(0, 0, 40, 2);
        return make([L(20, 7, "cu", { V: 100 }), L(20, 4, "fr4"), L(20, 1, "cu", { V: 0 })], { physics: "elec", grid: 1, gridReal: 10, depth: 40 });
      },
    },
    {
      name: "Câble coaxial", desc: "Âme à 100 V, isolant PTFE, blindage à 0 V : champ et capacité linéique",
      build() {
        const { circle, L, make } = exampleKit();
        circle(0, 0, 1.5); circle(0, 0, 5); circle(0, 0, 6);
        return make([L(0, 0, "cu", { V: 100 }), L(3, 0, "ptfe"), L(5.5, 0, "cu", { V: 0 })], { physics: "elec", grid: 0.5, gridReal: 1, depth: 1000, airMargin: 0.5 });
      },
    },
  ];
  const example = (i = 0) => EXAMPLES[i].build();
  const sure = (msg) => (S.edges.length ? ask({ msg, buttons: [["Annuler", null], ["Continuer", true, 1]] }) : Promise.resolve(true));
  $("btnNew").onclick = () => sure("Effacer le dessin en cours ?").then((ok) => { if (ok) { loadProject({ S: Object.assign(defaults(), { grid: S.grid, gridReal: S.gridReal, customMats: S.customMats }) }); fit(); } });
  $("btnExample").onclick = () =>
    ask({ msg: "Choisissez un exemple :", list: true, buttons: EXAMPLES.map((ex, i) => [ex.name + "\n" + ex.desc, i]).concat([["Annuler", null]]) })
      .then((i) => { if (i == null) return; return sure("Remplacer le dessin en cours par « " + EXAMPLES[i].name + " » ?").then((ok) => { if (ok) { loadProject(example(i)); fit(); } }); });
  $("btnSave").onclick = () => {
    const json = JSON.stringify({ app: "FieldTool", version: 1, S, view });
    ask({
      msg: "Projet au format texte. Téléchargez-le, ou copiez-le pour le coller dans Notes ou Fichiers.", text: json, readOnly: true,
      buttons: [
        ["Télécharger", () => {
          const a = document.createElement("a");
          a.href = URL.createObjectURL(new Blob([json], { type: "application/json" }));
          a.download = "circuit-" + new Date().toISOString().slice(0, 16).replace(/[:T]/g, "-") + ".json";
          a.click();
          setTimeout(() => URL.revokeObjectURL(a.href), 2000);
        }],
        ["Copier", (done, b) => {
          const sel = () => { $("askTxt").focus(); $("askTxt").select(); };
          try { navigator.clipboard.writeText(json).then(() => (b.textContent = "Copié"), sel); } catch (e) { sel(); }
        }],
        ["Fermer", null, 1],
      ],
    });
  };
  function loadText(txt) {
    try { const o = JSON.parse(txt); if (!o.S || !Array.isArray(o.S.edges)) throw new Error("format"); loadProject(o); fit(); }
    catch (err) { toast("Projet illisible."); }
  }
  $("btnOpen").onclick = () =>
    ask({
      msg: "Choisissez un fichier de projet, ou collez ci-dessous le texte d'un projet enregistré.", text: "",
      buttons: [["Annuler", null], ["Choisir un fichier", (done) => { done(null); $("fileOpen").click(); }], ["Charger le texte", true, 1]],
    }).then((v) => { if (v && v.trim()) loadText(v); });
  $("fileOpen").onchange = (e) => {
    const f = e.target.files[0];
    if (!f) return;
    const rd = new FileReader();
    rd.onload = () => loadText(rd.result);
    rd.readAsText(f);
    e.target.value = "";
  };
  window.addEventListener("keydown", (e) => {
    if (/INPUT|SELECT|TEXTAREA/.test(e.target.tagName) || $("dlg").open || $("ask").open) return;
    const k = e.key.toLowerCase();
    if ((e.ctrlKey || e.metaKey) && k === "z") { e.preventDefault(); (e.shiftKey ? $("btnRedo") : $("btnUndo")).click(); }
    else if ((e.ctrlKey || e.metaKey) && k === "y") { e.preventDefault(); $("btnRedo").click(); }
    else if (k === "delete" || k === "backspace") $("btnDel").click();
    else if (k === "enter") solve();
    else if (k === "g") setDrawMode(S.drawMode === "free" ? "grid" : "free");
    else { const t = { d: "draw", s: "edit", e: "erase", m: "mat", p: "probe" }[k]; if (t) setTool(t); }
  });

  // ---------- démarrage ----------
  initGL();
  new ResizeObserver(resize).observe(wrap);
  window.addEventListener("resize", resize);
  let stored = null;
  try { stored = JSON.parse(localStorage.getItem(KEY)); } catch (e) { stored = null; }
  if (stored && stored.S && Array.isArray(stored.S.edges)) loadProject(stored); else { syncInputs(); updateUI(); }
  if (window.innerWidth < 820) $("panel").classList.add("hidden");
  let palPref = null;
  try { palPref = localStorage.getItem("fieldtool.palette"); } catch (e) { palPref = null; }
  setPalette(palPref ? palPref === "1" : window.innerWidth >= 600);
  resize();
  window.FTapp = { get S() { return S; }, get mesh() { return mesh; }, get result() { return result; }, loadProject, example, examples: EXAMPLES, solve, fit, setTool };
})();
