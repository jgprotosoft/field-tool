// Géométrie : primitives (segment, arc), reconnaissance des tracés main levée, graphe planaire
(function (G) {
  "use strict";
  const FT = (G.FT = G.FT || {});

  // Arc défini par a, b et "bulge" = tan(angle/4), > 0 = sens trigonométrique de a vers b
  FT.arcParams = function (e) {
    const dx = e.bx - e.ax, dy = e.by - e.ay, d = Math.hypot(dx, dy), b = e.bulge;
    const nx = -dy / d, ny = dx / d;
    const k = ((d / 2) * (1 - b * b)) / (2 * b);
    const cx = (e.ax + e.bx) / 2 + nx * k, cy = (e.ay + e.by) / 2 + ny * k;
    return { cx, cy, r: Math.hypot(e.ax - cx, e.ay - cy), a0: Math.atan2(e.ay - cy, e.ax - cx), sweep: 4 * Math.atan(b) };
  };

  FT.edgePoints = function (e, stepDeg) {
    if (e.type !== "arc") return [[e.ax, e.ay], [e.bx, e.by]];
    const p = FT.arcParams(e);
    const n = Math.max(2, Math.ceil(Math.abs(p.sweep) / ((stepDeg || 6) * Math.PI / 180)));
    const out = [[e.ax, e.ay]];
    for (let i = 1; i < n; i++) {
      const a = p.a0 + (p.sweep * i) / n;
      out.push([p.cx + p.r * Math.cos(a), p.cy + p.r * Math.sin(a)]);
    }
    out.push([e.bx, e.by]);
    return out;
  };

  FT.edgeLength = function (e) {
    if (e.type !== "arc") return Math.hypot(e.bx - e.ax, e.by - e.ay);
    const p = FT.arcParams(e);
    return Math.abs(p.sweep) * p.r;
  };

  function ptSeg(x, y, ax, ay, bx, by) {
    const dx = bx - ax, dy = by - ay, l2 = dx * dx + dy * dy;
    let t = l2 > 0 ? ((x - ax) * dx + (y - ay) * dy) / l2 : 0;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    return Math.hypot(x - ax - t * dx, y - ay - t * dy);
  }
  FT.ptSeg = ptSeg;

  FT.distToEdge = function (e, x, y) {
    const p = FT.edgePoints(e, 6);
    let d = Infinity;
    for (let i = 0; i + 1 < p.length; i++) d = Math.min(d, ptSeg(x, y, p[i][0], p[i][1], p[i + 1][0], p[i + 1][1]));
    return d;
  };

  FT.modelBBox = function (edges) {
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const e of edges)
      for (const p of FT.edgePoints(e, 10)) {
        if (p[0] < x0) x0 = p[0];
        if (p[0] > x1) x1 = p[0];
        if (p[1] < y0) y0 = p[1];
        if (p[1] > y1) y1 = p[1];
      }
    return { x0, y0, x1, y1 };
  };

  // ---------- Reconnaissance d'un tracé main levée ----------
  function pathLength(p) {
    let l = 0;
    for (let i = 1; i < p.length; i++) l += Math.hypot(p[i].x - p[i - 1].x, p[i].y - p[i - 1].y);
    return l;
  }

  function resample(pts, ds) {
    const out = [pts[0]];
    let acc = 0;
    for (let i = 1; i < pts.length; i++) {
      let a = pts[i - 1];
      const b = pts[i];
      let d = Math.hypot(b.x - a.x, b.y - a.y);
      while (d > 0 && acc + d >= ds) {
        const t = (ds - acc) / d;
        a = { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t };
        out.push(a);
        d = Math.hypot(b.x - a.x, b.y - a.y);
        acc = 0;
      }
      acc += d;
    }
    const l = pts[pts.length - 1], o = out[out.length - 1];
    if (out.length < 2 || Math.hypot(l.x - o.x, l.y - o.y) > ds * 0.3) out.push(l);
    else out[out.length - 1] = l;
    return out;
  }

  function fitCircle(p, i0, i1) {
    const n = i1 - i0 + 1;
    let mx = 0, my = 0;
    for (let i = i0; i <= i1; i++) { mx += p[i].x; my += p[i].y; }
    mx /= n; my /= n;
    let suu = 0, suv = 0, svv = 0, suuu = 0, svvv = 0, suvv = 0, svuu = 0;
    for (let i = i0; i <= i1; i++) {
      const u = p[i].x - mx, v = p[i].y - my;
      suu += u * u; suv += u * v; svv += v * v;
      suuu += u * u * u; svvv += v * v * v; suvv += u * v * v; svuu += v * u * u;
    }
    const det = suu * svv - suv * suv;
    if (Math.abs(det) < 1e-12 * (suu + svv) * (suu + svv) + 1e-300) return null;
    const r1 = 0.5 * (suuu + suvv), r2 = 0.5 * (svvv + svuu);
    const uc = (r1 * svv - r2 * suv) / det, vc = (suu * r2 - suv * r1) / det;
    const r = Math.sqrt(uc * uc + vc * vc + (suu + svv) / n);
    let ss = 0;
    for (let i = i0; i <= i1; i++) {
      const d = Math.hypot(p[i].x - mx - uc, p[i].y - my - vc) - r;
      ss += d * d;
    }
    return { cx: mx + uc, cy: my + vc, r, rms: Math.sqrt(ss / n) };
  }

  function rdp(p, i0, i1, tol, out) {
    let dm = 0, im = -1;
    for (let i = i0 + 1; i < i1; i++) {
      const d = ptSeg(p[i].x, p[i].y, p[i0].x, p[i0].y, p[i1].x, p[i1].y);
      if (d > dm) { dm = d; im = i; }
    }
    if (im >= 0 && dm > tol) {
      rdp(p, i0, im, tol, out);
      out.push(im);
      rdp(p, im, i1, tol, out);
    }
  }

  // Point le plus proche sur une arête (ligne ou arc), avec sa distance.
  function nearestOnEdge(e, x, y) {
    if (e.type === "arc") {
      const a = FT.arcParams(e);
      let t = Math.atan2(y - a.cy, x - a.cx) - a.a0;
      const s = a.sweep;
      // ramener l'angle dans le balayage de l'arc
      if (s > 0) { while (t < 0) t += 2 * Math.PI; while (t > 2 * Math.PI) t -= 2 * Math.PI; }
      else { while (t > 0) t -= 2 * Math.PI; while (t < -2 * Math.PI) t += 2 * Math.PI; }
      if (Math.abs(t) > Math.abs(s)) {
        const da = Math.hypot(x - e.ax, y - e.ay), db = Math.hypot(x - e.bx, y - e.by);
        return da < db ? { x: e.ax, y: e.ay, d: da } : { x: e.bx, y: e.by, d: db };
      }
      const px = a.cx + a.r * Math.cos(a.a0 + t), py = a.cy + a.r * Math.sin(a.a0 + t);
      return { x: px, y: py, d: Math.hypot(x - px, y - py) };
    }
    const dx = e.bx - e.ax, dy = e.by - e.ay, L2 = dx * dx + dy * dy;
    const u = L2 ? Math.max(0, Math.min(1, ((x - e.ax) * dx + (y - e.ay) * dy) / L2)) : 0;
    const px = e.ax + u * dx, py = e.ay + u * dy;
    return { x: px, y: py, d: Math.hypot(x - px, y - py) };
  }

  // Accroche un point libre : d'abord aux sommets existants, puis aux arêtes, sinon le laisse en place.
  FT.snapFree = function (p, edges, r, skip) {
    let best = null, bd = r;
    for (const e of edges)
      for (const v of [[e.ax, e.ay], [e.bx, e.by]]) {
        if (skip && v[0] === skip.x && v[1] === skip.y) continue;
        const d = Math.hypot(v[0] - p.x, v[1] - p.y);
        if (d < bd) { bd = d; best = { x: v[0], y: v[1] }; }
      }
    if (best) return best;
    if (!skip) {
      for (const e of edges) { const q = nearestOnEdge(e, p.x, p.y); if (q.d < bd) { bd = q.d; best = q; } }
      if (best) return { x: +best.x.toFixed(4), y: +best.y.toFixed(4) };
    }
    return { x: +p.x.toFixed(3), y: +p.y.toFixed(3) };
  };

  // Mode « patate » : le tracé est gardé tel quel, simplifié en polyligne (tolérance tol en mm).
  // Les extrémités s'accrochent aux sommets ou arêtes à moins de snapR ; un tracé qui revient
  // près de son départ est fermé.
  FT.freehandStroke = function (raw, o) {
    if (!raw || raw.length < 2) return [];
    const tol = o.tol, snapR = o.snapR, edges = o.edges || [];
    const len = pathLength(raw);
    if (len < 4 * tol) return [];
    let p = resample(raw, Math.max(tol / 2, len / 2000));
    // lissage léger du tremblement, extrémités conservées
    for (let pass = 0; pass < 2; pass++) {
      const q = p.map((v) => ({ x: v.x, y: v.y }));
      for (let i = 1; i < p.length - 1; i++) { q[i].x = (p[i - 1].x + 2 * p[i].x + p[i + 1].x) / 4; q[i].y = (p[i - 1].y + 2 * p[i].y + p[i + 1].y) / 4; }
      p = q;
    }
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const v of p) { x0 = Math.min(x0, v.x); y0 = Math.min(y0, v.y); x1 = Math.max(x1, v.x); y1 = Math.max(y1, v.y); }
    const diag = Math.hypot(x1 - x0, y1 - y0);
    const gap = Math.hypot(p[p.length - 1].x - p[0].x, p[p.length - 1].y - p[0].y);
    const closed = gap < Math.max(snapR, 0.15 * diag) && len > 3 * gap && diag > 3 * tol;
    if (closed) p[p.length - 1] = { x: p[0].x, y: p[0].y };
    const idx = [0];
    rdp(p, 0, p.length - 1, tol, idx);
    idx.push(p.length - 1);
    idx.sort((a, b) => a - b);
    let pts = idx.map((i) => ({ x: p[i].x, y: p[i].y }));
    if (closed && pts.length < 4) {
      // contour fermé presque plat : on garde au moins un triangle
      const k = [0, Math.floor(p.length / 3), Math.floor((2 * p.length) / 3)];
      pts = k.map((i) => ({ x: p[i].x, y: p[i].y })).concat([{ x: p[0].x, y: p[0].y }]);
    }
    // accrochage des extrémités
    const a = FT.snapFree(pts[0], edges, snapR);
    pts[0] = a;
    if (closed) pts[pts.length - 1] = { x: a.x, y: a.y };
    else pts[pts.length - 1] = FT.snapFree(pts[pts.length - 1], edges, snapR);
    for (let i = 1; i < pts.length - 1; i++) pts[i] = { x: +pts[i].x.toFixed(3), y: +pts[i].y.toFixed(3) };
    // suppression des segments minuscules
    const out = [pts[0]];
    for (let i = 1; i < pts.length; i++) {
      const last = out[out.length - 1];
      if (Math.hypot(pts[i].x - last.x, pts[i].y - last.y) >= tol * 0.5 || i === pts.length - 1) out.push(pts[i]);
    }
    if (out.length > 2 && Math.hypot(out[out.length - 1].x - out[out.length - 2].x, out[out.length - 1].y - out[out.length - 2].y) < tol * 0.5) out.splice(out.length - 2, 1);
    const res = [];
    for (let i = 0; i + 1 < out.length; i++) {
      const u = out[i], v = out[i + 1];
      if (u.x !== v.x || u.y !== v.y) res.push({ type: "line", ax: u.x, ay: u.y, bx: v.x, by: v.y });
    }
    return res;
  };

  // raw : [{x,y}] en coordonnées monde ; grid : pas de grille. Retourne une liste d'arêtes accrochées.
  FT.recognizeStroke = function (raw, grid) {
    if (!raw || raw.length < 2) return [];
    const len = pathLength(raw);
    if (len < grid * 0.6) return [];
    const win = Math.max(grid * 0.5, len / 40);
    const W = 4, ds = win / W;
    const p = resample(raw, ds);
    const n = p.length;
    // lissage léger du tremblement (extrémités conservées)
    for (let pass = 0; pass < 2; pass++) {
      let px = p[0].x, py = p[0].y;
      for (let i = 1; i < n - 1; i++) {
        const x = (px + 2 * p[i].x + p[i + 1].x) / 4, y = (py + 2 * p[i].y + p[i + 1].y) / 4;
        px = p[i].x; py = p[i].y;
        p[i] = { x, y };
      }
    }
    const snap = (q) => ({ x: Math.round(q.x / grid) * grid, y: Math.round(q.y / grid) * grid });
    const closed = len > 3 * grid && Math.hypot(p[0].x - p[n - 1].x, p[0].y - p[n - 1].y) < Math.max(grid, 0.12 * len);

    // coins : maxima locaux de l'angle de braquage
    const turn = new Float64Array(n);
    for (let i = W; i < n - W; i++) {
      const ux = p[i].x - p[i - W].x, uy = p[i].y - p[i - W].y;
      const vx = p[i + W].x - p[i].x, vy = p[i + W].y - p[i].y;
      turn[i] = Math.abs(Math.atan2(ux * vy - uy * vx, ux * vx + uy * vy));
    }
    const corners = [];
    const TH = (45 * Math.PI) / 180;
    for (let i = W; i < n - W; i++) {
      if (turn[i] < TH) continue;
      let isMax = true;
      for (let j = Math.max(0, i - W); j <= Math.min(n - 1, i + W); j++) {
        if (j < i && turn[j] >= turn[i]) isMax = false;
        if (j > i && turn[j] > turn[i]) isMax = false;
      }
      if (isMax) corners.push(i);
    }

    // cercle complet
    if (closed) {
      const c = fitCircle(p, 0, n - 1);
      if (c && c.rms / c.r < (corners.length === 0 ? 0.1 : 0.03)) {
        const cc = snap({ x: c.cx, y: c.cy });
        const r = Math.max(grid, Math.round(c.r / grid) * grid);
        return [
          { type: "arc", ax: cc.x - r, ay: cc.y, bx: cc.x + r, by: cc.y, bulge: 1 },
          { type: "arc", ax: cc.x + r, ay: cc.y, bx: cc.x - r, by: cc.y, bulge: 1 },
        ];
      }
    }

    const idx = [0].concat(corners, [n - 1]);
    const V = idx.map((i) => snap(p[i]));
    if (closed) V[V.length - 1] = V[0];
    const out = [];
    const addLine = (a, b) => {
      if (a.x !== b.x || a.y !== b.y) out.push({ type: "line", ax: a.x, ay: a.y, bx: b.x, by: b.y });
    };
    for (let k = 0; k + 1 < idx.length; k++) {
      const i0 = idx[k], i1 = idx[k + 1], A = V[k], B = V[k + 1];
      if (A.x === B.x && A.y === B.y) {
        if (closed && idx.length === 2) {
          // contour fermé sans coin net et non circulaire : polygone simplifié
          const mids = [];
          rdp(p, i0, i1, Math.max(0.03 * len, grid * 0.5), mids);
          let prev = A;
          for (const m of mids) { const q = snap(p[m]); addLine(prev, q); prev = q; }
          addLine(prev, B);
        }
        continue;
      }
      const a = p[i0], b = p[i1];
      const d = Math.hypot(b.x - a.x, b.y - a.y);
      let dev = 0, sdev = 0;
      for (let i = i0 + 1; i < i1; i++) {
        const cr = ((b.x - a.x) * (p[i].y - a.y) - (b.y - a.y) * (p[i].x - a.x)) / (d || 1);
        if (Math.abs(cr) > dev) { dev = Math.abs(cr); sdev = cr; }
      }
      if (dev < 0.1 * d) { addLine(A, B); continue; }
      const c = fitCircle(p, i0, i1);
      if (c && c.rms / c.r < 0.15) {
        let deg = (4 * Math.atan((2 * dev) / d) * 180) / Math.PI;
        if (Math.abs(deg - 90) < 25) deg = 90;
        else if (Math.abs(deg - 180) < 30) deg = 180;
        else deg = Math.min(330, Math.max(30, Math.round(deg / 15) * 15));
        const bulge = Math.tan((deg * Math.PI) / 720) * (sdev < 0 ? 1 : -1);
        out.push({ type: "arc", ax: A.x, ay: A.y, bx: B.x, by: B.y, bulge });
        continue;
      }
      // forme quelconque : polyligne simplifiée
      const mids = [];
      rdp(p, i0, i1, Math.max(0.08 * d, grid * 0.5), mids);
      let prev = A;
      for (const m of mids) { const q = snap(p[m]); addLine(prev, q); prev = q; }
      addLine(prev, B);
    }
    return out;
  };

  // ---------- Graphe planaire : découpe aux intersections, fusion des points ----------
  // box : {x0,y0,x1,y1} boîte d'air ajoutée comme frontière extérieure
  FT.buildPSLG = function (edges, box) {
    const S = [];
    edges.forEach((e, idx) => {
      const p = FT.edgePoints(e, 6);
      for (let i = 0; i + 1 < p.length; i++)
        if (p[i][0] !== p[i + 1][0] || p[i][1] !== p[i + 1][1]) S.push([p[i][0], p[i][1], p[i + 1][0], p[i + 1][1], idx, 0]);
    });
    S.push([box.x0, box.y0, box.x1, box.y0, -1, 1], [box.x1, box.y0, box.x1, box.y1, -2, 1],
      [box.x1, box.y1, box.x0, box.y1, -3, 1], [box.x0, box.y1, box.x0, box.y0, -4, 1]);
    const diag = Math.hypot(box.x1 - box.x0, box.y1 - box.y0);
    const tol = 1e-7 * diag;
    const ns = S.length;
    const ts = S.map(() => [0, 1]);
    const onSeg = (px, py, s, list) => {
      const rx = s[2] - s[0], ry = s[3] - s[1], l2 = rx * rx + ry * ry;
      const t = ((px - s[0]) * rx + (py - s[1]) * ry) / l2;
      const te = tol / Math.sqrt(l2);
      if (t <= te || t >= 1 - te) return;
      if (Math.hypot(px - s[0] - t * rx, py - s[1] - t * ry) < tol) list.push(t);
    };
    for (let i = 0; i < ns; i++) {
      const a = S[i];
      const ax0 = Math.min(a[0], a[2]) - tol, ax1 = Math.max(a[0], a[2]) + tol;
      const ay0 = Math.min(a[1], a[3]) - tol, ay1 = Math.max(a[1], a[3]) + tol;
      for (let j = i + 1; j < ns; j++) {
        const b = S[j];
        if (Math.min(b[0], b[2]) > ax1 || Math.max(b[0], b[2]) < ax0 || Math.min(b[1], b[3]) > ay1 || Math.max(b[1], b[3]) < ay0) continue;
        onSeg(b[0], b[1], a, ts[i]); onSeg(b[2], b[3], a, ts[i]);
        onSeg(a[0], a[1], b, ts[j]); onSeg(a[2], a[3], b, ts[j]);
        const rx = a[2] - a[0], ry = a[3] - a[1], sx = b[2] - b[0], sy = b[3] - b[1];
        const den = rx * sy - ry * sx;
        if (Math.abs(den) > 1e-12 * Math.hypot(rx, ry) * Math.hypot(sx, sy)) {
          const qx = b[0] - a[0], qy = b[1] - a[1];
          const t = (qx * sy - qy * sx) / den, u = (qx * ry - qy * rx) / den;
          if (t > 0 && t < 1 && u > 0 && u < 1) { ts[i].push(t); ts[j].push(u); }
        }
      }
    }
    // registre de nœuds avec fusion
    const px = [], py = [];
    const mt = 2 * tol, cs = 8 * mt;
    const hash = new Map();
    const nodeId = (x, y) => {
      const ix = Math.floor(x / cs), iy = Math.floor(y / cs);
      for (let i = ix - 1; i <= ix + 1; i++)
        for (let j = iy - 1; j <= iy + 1; j++) {
          const l = hash.get(i + "," + j);
          if (l) for (const id of l) if (Math.abs(px[id] - x) <= mt && Math.abs(py[id] - y) <= mt) return id;
        }
      const id = px.length;
      px.push(x); py.push(y);
      const key = ix + "," + iy;
      const l = hash.get(key);
      if (l) l.push(id); else hash.set(key, [id]);
      return id;
    };
    const segMap = new Map();
    for (let i = 0; i < ns; i++) {
      const s = S[i];
      const t = ts[i].sort((a, b) => a - b);
      let prev = nodeId(s[0], s[1]);
      for (let k = 1; k < t.length; k++) {
        const id = t[k] === 1 ? nodeId(s[2], s[3]) : nodeId(s[0] + (s[2] - s[0]) * t[k], s[1] + (s[3] - s[1]) * t[k]);
        if (id !== prev) {
          const key = prev < id ? prev + "_" + id : id + "_" + prev;
          const ex = segMap.get(key);
          if (!ex) segMap.set(key, { a: prev, b: id, parent: s[4], outer: s[5] });
          else if (s[5]) ex.outer = 1;
        }
        prev = id;
      }
    }
    return { px, py, segs: Array.from(segMap.values()), box };
  };
})(typeof window !== "undefined" ? window : globalThis);
