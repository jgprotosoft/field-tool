// Mailleur triangulaire : Delaunay incrémental conforme, taille adaptée aux entrefers
(function (G) {
  "use strict";
  const FT = (G.FT = G.FT || {});

  // ---------- Delaunay incrémental (Bowyer-Watson) ----------
  function Del(x0, y0, x1, y1) {
    this.px = [x0, x1, x1, x0];
    this.py = [y0, y0, y1, y1];
    this.tv = [0, 1, 2, 0, 2, 3];
    this.tn = [-1, 1, -1, -1, -1, 0]; // voisin en face du sommet k
    this.alive = [1, 1];
    this.mark = [0, 0];
    this.free = [];
    this.stamp = 0;
    this.last = 0;
    this.sm = [];
    this.em = [];
  }
  Del.prototype.orient = function (a, b, x, y) {
    const px = this.px, py = this.py;
    return (px[b] - px[a]) * (y - py[a]) - (py[b] - py[a]) * (x - px[a]);
  };
  Del.prototype.inCircle = function (t, x, y) {
    const px = this.px, py = this.py, tv = this.tv;
    const a = tv[3 * t], b = tv[3 * t + 1], c = tv[3 * t + 2];
    const ax = px[a] - x, ay = py[a] - y, bx = px[b] - x, by = py[b] - y, cx = px[c] - x, cy = py[c] - y;
    const a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const det = a2 * (bx * cy - cx * by) - b2 * (ax * cy - cx * ay) + c2 * (ax * by - bx * ay);
    return det > 1e-12 * (a2 * b2 + b2 * c2 + c2 * a2);
  };
  Del.prototype.locate = function (x, y) {
    const tv = this.tv, tn = this.tn;
    let t = this.last;
    const lim = this.alive.length * 2 + 200;
    let it = 0;
    outer: while (it++ < lim) {
      const k0 = it % 3;
      for (let kk = 0; kk < 3; kk++) {
        const k = (k0 + kk) % 3;
        const a = tv[3 * t + ((k + 1) % 3)], b = tv[3 * t + ((k + 2) % 3)];
        if (this.orient(a, b, x, y) < -1e-14) {
          const n = tn[3 * t + k];
          if (n < 0) return -1;
          t = n;
          continue outer;
        }
      }
      return t;
    }
    for (t = 0; t < this.alive.length; t++) {
      if (!this.alive[t]) continue;
      let ok = true;
      for (let k = 0; k < 3 && ok; k++)
        if (this.orient(tv[3 * t + ((k + 1) % 3)], tv[3 * t + ((k + 2) % 3)], x, y) < -1e-14) ok = false;
      if (ok) return t;
    }
    return -1;
  };
  Del.prototype.insert = function (x, y) {
    const px = this.px, py = this.py, tv = this.tv, tn = this.tn, mark = this.mark;
    const p = px.length;
    px.push(x); py.push(y);
    const t = this.locate(x, y);
    if (t < 0) return -1;
    const st = ++this.stamp;
    const cav = [t];
    mark[t] = st;
    for (let q = 0; q < cav.length; q++) {
      const c = cav[q];
      for (let k = 0; k < 3; k++) {
        const n = tn[3 * c + k];
        if (n >= 0 && mark[n] !== st && this.inCircle(n, x, y)) { mark[n] = st; cav.push(n); }
      }
    }
    let bnd, guard = 0;
    for (;;) {
      bnd = [];
      let bad = -1, badN = -1;
      scan: for (let q = 0; q < cav.length; q++) {
        const c = cav[q];
        if (mark[c] !== st) continue;
        for (let k = 0; k < 3; k++) {
          const n = tn[3 * c + k];
          if (n < 0 || mark[n] !== st) {
            const a = tv[3 * c + ((k + 1) % 3)], b = tv[3 * c + ((k + 2) % 3)];
            const o = this.orient(a, b, x, y);
            if (o <= 1e-10 * Math.hypot(px[b] - px[a], py[b] - py[a])) { bad = c; badN = n; break scan; }
            bnd.push(a, b, n, c);
          }
        }
      }
      if (bad < 0) break;
      if (++guard > 60) return -1;
      if (bad === t) {
        if (badN < 0) return -1;
        mark[badN] = st;
        cav.push(badN);
      } else mark[bad] = -1;
    }
    const news = [];
    for (let i = 0; i < bnd.length; i += 4) {
      const a = bnd[i], b = bnd[i + 1], n = bnd[i + 2], c = bnd[i + 3];
      let T;
      if (this.free.length) T = this.free.pop();
      else { T = this.alive.length; this.alive.push(1); mark.push(0); tv.push(0, 0, 0); tn.push(-1, -1, -1); }
      tv[3 * T] = a; tv[3 * T + 1] = b; tv[3 * T + 2] = p;
      tn[3 * T + 2] = n;
      this.alive[T] = 1; mark[T] = 0;
      if (n >= 0) for (let k = 0; k < 3; k++) if (tn[3 * n + k] === c) { tn[3 * n + k] = T; break; }
      this.sm[a] = T; this.em[b] = T;
      news.push(T);
    }
    for (const T of news) { tn[3 * T] = this.sm[tv[3 * T + 1]]; tn[3 * T + 1] = this.em[tv[3 * T]]; }
    for (const c of cav) if (mark[c] === st) { this.alive[c] = 0; this.free.push(c); }
    this.last = news[0];
    return p;
  };
  FT._Del = Del;

  // ---------- Génération du maillage ----------
  // pslg : sortie de buildPSLG ; o : {hmax (mm), refine (hmax/hmin), grad}
  FT.generateMesh = function (pslg, o) {
    const T0 = Date.now();
    const bx = pslg.box;
    const sc = 1 / Math.max(bx.x1 - bx.x0, bx.y1 - bx.y0);
    const X = pslg.px.map((v) => (v - bx.x0) * sc), Y = pslg.py.map((v) => (v - bx.y0) * sc);
    const BW = (bx.x1 - bx.x0) * sc, BH = (bx.y1 - bx.y0) * sc;
    const hmax = o.hmax * sc, hmin = hmax / (o.refine || 32), hfar = hmax * 4, g = o.grad || 0.3;
    const warnings = [];

    // 1. découpage initial des segments
    const pa = [], pb = [], ppar = [], pout = [];
    for (const s of pslg.segs) {
      const L = Math.hypot(X[s.b] - X[s.a], Y[s.b] - Y[s.a]);
      const n = Math.max(1, Math.ceil(L / (s.outer ? hfar : hmax) - 1e-9));
      let prev = s.a;
      for (let k = 1; k <= n; k++) {
        let id = s.b;
        if (k < n) { id = X.length; X.push(X[s.a] + ((X[s.b] - X[s.a]) * k) / n); Y.push(Y[s.a] + ((Y[s.b] - Y[s.a]) * k) / n); }
        pa.push(prev); pb.push(id); ppar.push(s.parent); pout.push(s.outer ? 1 : 0);
        prev = id;
      }
    }
    const np = pa.length;
    // 2. taille locale d'après la distance aux autres frontières (entrefers, parois minces)
    const mx = new Float64Array(np), my = new Float64Array(np), hl = new Float64Array(np), ss = new Float64Array(np);
    for (let i = 0; i < np; i++) {
      mx[i] = (X[pa[i]] + X[pb[i]]) / 2; my[i] = (Y[pa[i]] + Y[pb[i]]) / 2;
      hl[i] = Math.hypot(X[pb[i]] - X[pa[i]], Y[pb[i]] - Y[pa[i]]) / 2;
      ss[i] = pout[i] ? hfar : hmax;
      // taille imposée par l'appelant (épaisseur de peau des conducteurs), en mm
      if (o.cap && !pout[i]) {
        const c = o.cap(X[pa[i]] / sc + bx.x0, Y[pa[i]] / sc + bx.y0, X[pb[i]] / sc + bx.x0, Y[pb[i]] / sc + bx.y0) * sc;
        if (c < ss[i]) ss[i] = Math.max(hmin, c);
      }
    }
    const ptSeg = FT.ptSeg;
    const lim = hmax / 0.6;
    for (let i = 0; i < np; i++) {
      if (pout[i]) continue;
      for (let j = i + 1; j < np; j++) {
        if (pout[j] || ppar[i] === ppar[j]) continue;
        const dm = Math.hypot(mx[i] - mx[j], my[i] - my[j]) - hl[i] - hl[j];
        if (dm > lim) continue;
        const a = pa[i], b = pb[i], c = pa[j], d = pb[j];
        if (a === c || a === d || b === c || b === d) continue;
        const dist = Math.min(ptSeg(X[a], Y[a], X[c], Y[c], X[d], Y[d]), ptSeg(X[b], Y[b], X[c], Y[c], X[d], Y[d]),
          ptSeg(X[c], Y[c], X[a], Y[a], X[b], Y[b]), ptSeg(X[d], Y[d], X[a], Y[a], X[b], Y[b]));
        const s = Math.max(hmin, 0.6 * dist);
        if (s < ss[i]) ss[i] = s;
        if (s < ss[j]) ss[j] = s;
      }
    }
    // 3. gradation
    const order = Array.from({ length: np }, (_, i) => i).sort((i, j) => ss[i] - ss[j]);
    for (const i of order) {
      const si = ss[i];
      if (si >= hmax) break;
      for (let j = 0; j < np; j++) {
        const v = si + g * Math.hypot(mx[i] - mx[j], my[i] - my[j]);
        if (v < ss[j]) ss[j] = v;
      }
    }
    const hf = (x, y) => {
      let h = hfar;
      for (let i = 0; i < np; i++) {
        const v = ss[i] + g * Math.hypot(x - mx[i], y - my[i]);
        if (v < h) h = v;
      }
      return h;
    };
    // 4. subdivision finale des frontières
    let fa = [], fb = [], fo = [];
    for (let i = 0; i < np; i++) {
      const n = Math.max(1, Math.round((2 * hl[i]) / ss[i] + 0.3));
      let prev = pa[i];
      for (let k = 1; k <= n; k++) {
        let id = pb[i];
        if (k < n) { id = X.length; X.push(X[pa[i]] + ((X[pb[i]] - X[pa[i]]) * k) / n); Y.push(Y[pa[i]] + ((Y[pb[i]] - Y[pa[i]]) * k) / n); }
        fa.push(prev); fb.push(id); fo.push(pout[i]);
        prev = id;
      }
    }
    // 5. triangulation : nœuds frontière puis nœuds intérieurs (réseau hexagonal hiérarchique)
    const mrg = 0.3;
    const D = new Del(-mrg, -mrg, BW + mrg, BH + mrg);
    const map = new Int32Array(X.length); // indice pslg -> indice Delaunay
    for (let i = 0; i < X.length; i++) {
      let id = D.insert(X[i], Y[i]);
      if (id < 0) id = D.insert(X[i] + 1e-9, Y[i] + 1.3e-9);
      if (id < 0) warnings.push("nœud frontière non inséré");
      map[i] = id;
    }
    for (let i = 0; i < fa.length; i++) { fa[i] = map[fa[i]]; fb[i] = map[fb[i]]; }
    // grille de recherche des segments frontière
    const cs = hfar, gnx = Math.ceil(BW / cs) + 1, gny = Math.ceil(BH / cs) + 1;
    const cells = Array.from({ length: gnx * gny }, () => []);
    const cellOf = (v) => Math.max(0, Math.floor(v / cs));
    for (let i = 0; i < fa.length; i++) {
      if (fa[i] < 0 || fb[i] < 0) continue;
      const x0 = D.px[fa[i]], y0 = D.py[fa[i]], x1 = D.px[fb[i]], y1 = D.py[fb[i]];
      const i0 = cellOf(Math.min(x0, x1)), i1 = Math.min(gnx - 1, cellOf(Math.max(x0, x1)));
      const j0 = cellOf(Math.min(y0, y1)), j1 = Math.min(gny - 1, cellOf(Math.max(y0, y1)));
      for (let ii = i0; ii <= i1; ii++) for (let jj = j0; jj <= j1; jj++) cells[jj * gnx + ii].push(i);
    }
    const nearB = (x, y, r) => {
      const i0 = cellOf(x - r), i1 = Math.min(gnx - 1, cellOf(x + r)), j0 = cellOf(y - r), j1 = Math.min(gny - 1, cellOf(y + r));
      for (let ii = i0; ii <= i1; ii++)
        for (let jj = j0; jj <= j1; jj++)
          for (const s of cells[jj * gnx + ii])
            if (ptSeg(x, y, D.px[fa[s]], D.py[fa[s]], D.px[fb[s]], D.py[fb[s]]) < r) return true;
      return false;
    };
    const R3 = Math.sqrt(3) / 2;
    const maxPts = o.maxNodes || 150000;
    let nInt = 0;
    const visit = (x, y, L) => {
      const s = hmin * (1 << L);
      if (x < -s || x > BW + s || y < -s || y > BH + s) return;
      const h = hf(x, y);
      if (L > 0 && h < 0.75 * s) {
        const e = s / 2;
        visit(x, y, L - 1); visit(x + e, y, L - 1); visit(x + e / 2, y + e * R3, L - 1); visit(x - e / 2, y + e * R3, L - 1);
        return;
      }
      if (x <= 0 || x >= BW || y <= 0 || y >= BH || nInt > maxPts) return;
      if (nearB(x, y, 0.5 * Math.min(h, 1.5 * s))) return;
      if (D.insert(x, y) >= 0) nInt++;
    };
    let Ltop = 0;
    while (hmin * (1 << Ltop) < hfar * 0.99) Ltop++;
    const s0 = hmin * (1 << Ltop);
    const off = 0.137 * s0; // décalage pour éviter les alignements avec la grille de dessin
    for (let j = -1; j * s0 * R3 <= BH + s0; j++)
      for (let i = -1 - Math.ceil(j / 2); i * s0 + (j * s0) / 2 <= BW + s0; i++) visit(i * s0 + (j * s0) / 2 + off, j * s0 * R3 + off * 0.61, Ltop);
    if (nInt > maxPts) warnings.push("nombre de nœuds limité : réduire la finesse");

    // 6. conformité : découper les segments absents de la triangulation
    const NN = 1 << 26;
    const key = (a, b) => (a < b ? a * NN + b : b * NN + a);
    let missing = 0;
    for (let iter = 0; iter < 30; iter++) {
      const es = new Set();
      const tv = D.tv;
      for (let t = 0; t < D.alive.length; t++) {
        if (!D.alive[t]) continue;
        es.add(key(tv[3 * t], tv[3 * t + 1])); es.add(key(tv[3 * t + 1], tv[3 * t + 2])); es.add(key(tv[3 * t + 2], tv[3 * t]));
      }
      const na = [], nb = [], no = [];
      missing = 0;
      for (let i = 0; i < fa.length; i++) {
        const a = fa[i], b = fb[i];
        if (a < 0 || b < 0) continue;
        if (es.has(key(a, b))) { na.push(a); nb.push(b); no.push(fo[i]); continue; }
        missing++;
        const m = D.insert((D.px[a] + D.px[b]) / 2, (D.py[a] + D.py[b]) / 2);
        if (m < 0) { na.push(a); nb.push(b); no.push(fo[i]); continue; }
        na.push(a, m); nb.push(m, b); no.push(fo[i], fo[i]);
      }
      fa = na; fb = nb; fo = no;
      if (!missing) break;
    }
    if (missing) warnings.push(missing + " segment(s) de frontière non respecté(s)");

    // 7. régions par remplissage, extraction
    const cons = new Set();
    for (let i = 0; i < fa.length; i++) cons.add(key(fa[i], fb[i]));
    const nT = D.alive.length, tv = D.tv, tn = D.tn;
    const reg = new Int32Array(nT).fill(-1);
    let nReg = 0;
    const outside = [];
    for (let t = 0; t < nT; t++) {
      if (!D.alive[t] || reg[t] >= 0) continue;
      const r = nReg++;
      let isOut = false;
      const stack = [t];
      reg[t] = r;
      while (stack.length) {
        const c = stack.pop();
        if (tv[3 * c] < 4 || tv[3 * c + 1] < 4 || tv[3 * c + 2] < 4) isOut = true;
        for (let k = 0; k < 3; k++) {
          const n = tn[3 * c + k];
          if (n < 0 || reg[n] >= 0) continue;
          if (cons.has(key(tv[3 * c + ((k + 1) % 3)], tv[3 * c + ((k + 2) % 3)]))) continue;
          reg[n] = r;
          stack.push(n);
        }
      }
      outside.push(isOut);
    }
    const regMap = new Int32Array(nReg).fill(-1);
    let nr = 0;
    for (let r = 0; r < nReg; r++) if (!outside[r]) regMap[r] = nr++;
    const nodeMap = new Int32Array(D.px.length).fill(-1);
    const ox = [], oy = [], otri = [], oreg = [];
    const nid = (v) => {
      if (nodeMap[v] < 0) { nodeMap[v] = ox.length; ox.push(D.px[v] / sc + bx.x0); oy.push(D.py[v] / sc + bx.y0); }
      return nodeMap[v];
    };
    for (let t = 0; t < nT; t++) {
      if (!D.alive[t] || regMap[reg[t]] < 0) continue;
      otri.push(nid(tv[3 * t]), nid(tv[3 * t + 1]), nid(tv[3 * t + 2]));
      oreg.push(regMap[reg[t]]);
    }
    const outer = new Uint8Array(ox.length);
    const bseg = [];
    for (let i = 0; i < fa.length; i++) {
      const a = nodeMap[fa[i]], b = nodeMap[fb[i]];
      if (a < 0 || b < 0) continue;
      if (fo[i]) { outer[a] = 1; outer[b] = 1; }
      else bseg.push(a, b);
    }
    return {
      x: Float64Array.from(ox), y: Float64Array.from(oy), tri: Int32Array.from(otri), reg: Int32Array.from(oreg),
      nReg: nr, outer, bseg: Int32Array.from(bseg), box: bx, warnings, ms: Date.now() - T0,
    };
  };

  // ---------- Localisation d'un point dans le maillage ----------
  FT.buildLocator = function (mesh) {
    const nt = mesh.tri.length / 3, bx = mesh.box;
    const n = Math.max(4, Math.ceil(Math.sqrt(nt / 2)));
    const w = (bx.x1 - bx.x0) / n, h = (bx.y1 - bx.y0) / n;
    const cells = Array.from({ length: n * n }, () => []);
    const cl = (v) => (v < 0 ? 0 : v >= n ? n - 1 : v);
    const { x, y, tri } = mesh;
    for (let t = 0; t < nt; t++) {
      const a = tri[3 * t], b = tri[3 * t + 1], c = tri[3 * t + 2];
      const i0 = cl(Math.floor((Math.min(x[a], x[b], x[c]) - bx.x0) / w)), i1 = cl(Math.floor((Math.max(x[a], x[b], x[c]) - bx.x0) / w));
      const j0 = cl(Math.floor((Math.min(y[a], y[b], y[c]) - bx.y0) / h)), j1 = cl(Math.floor((Math.max(y[a], y[b], y[c]) - bx.y0) / h));
      for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) cells[j * n + i].push(t);
    }
    return {
      // retourne l'indice du triangle contenant (px,py) ou -1 ; bary reçoit les coordonnées barycentriques
      find(px, py, bary) {
        if (px < bx.x0 || px > bx.x1 || py < bx.y0 || py > bx.y1) return -1;
        const list = cells[cl(Math.floor((py - bx.y0) / h)) * n + cl(Math.floor((px - bx.x0) / w))];
        for (const t of list) {
          const a = tri[3 * t], b = tri[3 * t + 1], c = tri[3 * t + 2];
          const d = (x[b] - x[a]) * (y[c] - y[a]) - (x[c] - x[a]) * (y[b] - y[a]);
          const l1 = ((x[b] - px) * (y[c] - py) - (x[c] - px) * (y[b] - py)) / d;
          const l2 = ((x[c] - px) * (y[a] - py) - (x[a] - px) * (y[c] - py)) / d;
          const l3 = 1 - l1 - l2;
          if (l1 >= -1e-9 && l2 >= -1e-9 && l3 >= -1e-9) {
            if (bary) { bary[0] = l1; bary[1] = l2; bary[2] = l3; }
            return t;
          }
        }
        return -1;
      },
    };
  };
})(typeof window !== "undefined" ? window : globalThis);
