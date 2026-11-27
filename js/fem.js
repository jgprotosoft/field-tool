// Solveur éléments finis 2D plan (triangles linéaires)
//  - magnétostatique : ∇·(ν∇A) = −J − rot(ν·Br), Newton-Raphson pour les matériaux saturables
//  - magnétodynamique sinusoïdale : ∇·(ν∇A) − jωσA = −J (courants induits, effet de peau)
//  - électrostatique : ∇·(ε∇V) = −ρ
(function (G) {
  "use strict";
  const FT = (G.FT = G.FT || {});
  const MU0 = FT.MU0, EPS0 = FT.EPS0;

  // Géométrie élémentaire en mètres : gradients des fonctions de forme (b,c)/(2Δ)
  function elemGeom(mesh) {
    const nt = mesh.tri.length / 3, { x, y, tri } = mesh;
    const gb = new Float64Array(3 * nt), gc = new Float64Array(3 * nt), area = new Float64Array(nt);
    for (let e = 0; e < nt; e++) {
      const i = tri[3 * e], j = tri[3 * e + 1], k = tri[3 * e + 2];
      const xi = x[i] * 1e-3, yi = y[i] * 1e-3, xj = x[j] * 1e-3, yj = y[j] * 1e-3, xk = x[k] * 1e-3, yk = y[k] * 1e-3;
      gb[3 * e] = yj - yk; gb[3 * e + 1] = yk - yi; gb[3 * e + 2] = yi - yj;
      gc[3 * e] = xk - xj; gc[3 * e + 1] = xi - xk; gc[3 * e + 2] = xj - xi;
      area[e] = 0.5 * ((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi));
    }
    return { gb, gc, area };
  }
  FT.elemGeom = elemGeom;

  // Structure creuse CSR pour les degrés de liberté libres
  function buildStructure(tri, dof, nd) {
    const nt = tri.length / 3;
    const rows = Array.from({ length: nd }, () => []);
    for (let e = 0; e < nt; e++)
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        const r = rows[i];
        for (let b = 0; b < 3; b++) {
          const j = dof[tri[3 * e + b]];
          if (j >= 0 && r.indexOf(j) < 0) r.push(j);
        }
      }
    const rp = new Int32Array(nd + 1);
    for (let i = 0; i < nd; i++) rp[i + 1] = rp[i] + rows[i].length;
    const ci = new Int32Array(rp[nd]), dg = new Int32Array(nd);
    for (let i = 0; i < nd; i++) {
      const r = rows[i].sort((p, q) => p - q);
      for (let k = 0; k < r.length; k++) { ci[rp[i] + k] = r[k]; if (r[k] === i) dg[i] = rp[i] + k; }
    }
    const pos = new Int32Array(9 * nt).fill(-1);
    for (let e = 0; e < nt; e++)
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        for (let b = 0; b < 3; b++) {
          const j = dof[tri[3 * e + b]];
          if (j < 0) continue;
          let lo = rp[i], hi = rp[i + 1] - 1;
          while (lo < hi) { const m = (lo + hi) >> 1; if (ci[m] < j) lo = m + 1; else hi = m; }
          pos[9 * e + 3 * a + b] = lo;
        }
      }
    return { rp, ci, dg, pos, va: new Float64Array(rp[nd]), nd };
  }

  // Gradient conjugué préconditionné Gauss-Seidel symétrique
  function pcg(M, b, x, tol, maxit) {
    const { rp, ci, va, dg, nd: n } = M;
    const r = new Float64Array(n), z = new Float64Array(n), p = new Float64Array(n), q = new Float64Array(n);
    const mul = (v, out) => {
      for (let i = 0; i < n; i++) { let s = 0; for (let k = rp[i]; k < rp[i + 1]; k++) s += va[k] * v[ci[k]]; out[i] = s; }
    };
    const prec = () => {
      for (let i = 0; i < n; i++) {
        let s = r[i];
        for (let k = rp[i]; k < dg[i]; k++) s -= va[k] * z[ci[k]];
        z[i] = s / va[dg[i]];
      }
      for (let i = n - 1; i >= 0; i--) {
        let s = 0;
        for (let k = dg[i] + 1; k < rp[i + 1]; k++) s += va[k] * z[ci[k]];
        z[i] -= s / va[dg[i]];
      }
    };
    mul(x, q);
    let bn = 0, rn = 0;
    for (let i = 0; i < n; i++) { r[i] = b[i] - q[i]; bn += b[i] * b[i]; rn += r[i] * r[i]; }
    if (bn === 0) { x.fill(0); return { it: 0, res: 0 }; }
    const target = tol * tol * bn;
    if (rn <= target) return { it: 0, res: Math.sqrt(rn / bn) };
    prec();
    p.set(z);
    let rz = 0;
    for (let i = 0; i < n; i++) rz += r[i] * z[i];
    let it = 0;
    for (; it < maxit; it++) {
      mul(p, q);
      let pq = 0;
      for (let i = 0; i < n; i++) pq += p[i] * q[i];
      const al = rz / pq;
      rn = 0;
      for (let i = 0; i < n; i++) { x[i] += al * p[i]; r[i] -= al * q[i]; rn += r[i] * r[i]; }
      if (rn <= target) { it++; break; }
      prec();
      let rz2 = 0;
      for (let i = 0; i < n; i++) rz2 += r[i] * z[i];
      const be = rz2 / rz;
      rz = rz2;
      for (let i = 0; i < n; i++) p[i] = z[i] + be * p[i];
    }
    return { it, res: Math.sqrt(rn / bn) };
  }

  // props[r] = {mur, Js, J (A/m²), brx, bry (T)} ; frontière extérieure A = 0
  FT.solveMagnetic = function (mesh, props) {
    const T0 = Date.now();
    const n = mesh.x.length, tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg;
    const { gb, gc, area } = elemGeom(mesh);
    const dof = new Int32Array(n);
    let nd = 0;
    for (let i = 0; i < n; i++) dof[i] = mesh.outer[i] ? -1 : nd++;
    const M = buildStructure(tri, dof, nd);
    const nuE = new Float64Array(nt), dnuE = new Float64Array(nt);
    const f = new Float64Array(nd);
    let nonlinear = false;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      nuE[e] = 1 / (MU0 * p.mur);
      if (p.Js > 0) nonlinear = true;
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        f[i] += (p.J * area[e]) / 3 + (nuE[e] * (p.brx * gc[3 * e + a] - p.bry * gb[3 * e + a])) / 2;
      }
    }
    const A = new Float64Array(n);
    const gxE = new Float64Array(nt), gyE = new Float64Array(nt);
    const update = (Av) => {
      for (let e = 0; e < nt; e++) {
        const a0 = Av[tri[3 * e]], a1 = Av[tri[3 * e + 1]], a2 = Av[tri[3 * e + 2]], d = 2 * area[e];
        const gx = (a0 * gb[3 * e] + a1 * gb[3 * e + 1] + a2 * gb[3 * e + 2]) / d;
        const gy = (a0 * gc[3 * e] + a1 * gc[3 * e + 1] + a2 * gc[3 * e + 2]) / d;
        gxE[e] = gx; gyE[e] = gy;
        const p = props[reg[e]];
        if (p.Js > 0) { const r = FT.nuOfB(p.mur, p.Js, Math.hypot(gx, gy)); nuE[e] = r[0]; dnuE[e] = r[1]; }
      }
    };
    const residual = (out) => {
      out.set(f);
      for (let e = 0; e < nt; e++) {
        const nu = nuE[e], gx = gxE[e], gy = gyE[e];
        for (let a = 0; a < 3; a++) {
          const i = dof[tri[3 * e + a]];
          if (i >= 0) out[i] -= (nu * (gx * gb[3 * e + a] + gy * gc[3 * e + a])) / 2;
        }
      }
      let s = 0;
      for (let i = 0; i < nd; i++) s += out[i] * out[i];
      return Math.sqrt(s);
    };
    const assemble = () => {
      const va = M.va, pos = M.pos;
      va.fill(0);
      for (let e = 0; e < nt; e++) {
        const nu = nuE[e], dn = dnuE[e], ar = area[e], gx = gxE[e], gy = gyE[e];
        for (let a = 0; a < 3; a++) {
          const qa = (gx * gb[3 * e + a] + gy * gc[3 * e + a]) / 2;
          for (let b = 0; b < 3; b++) {
            const k = pos[9 * e + 3 * a + b];
            if (k < 0) continue;
            let v = (nu * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
            if (dn > 0) v += ((2 * dn) / ar) * qa * ((gx * gb[3 * e + b] + gy * gc[3 * e + b]) / 2);
            va[k] += v;
          }
        }
      }
    };
    let fn = 0;
    for (let i = 0; i < nd; i++) fn += f[i] * f[i];
    fn = Math.sqrt(fn);
    const r = new Float64Array(nd), dx = new Float64Array(nd), At = new Float64Array(n);
    let iters = 0, cgIt = 0, converged = true, relRes = 0;
    if (fn > 0) {
      update(A);
      let rn = residual(r);
      const maxNewton = nonlinear ? 40 : 1;
      converged = !nonlinear;
      for (let it = 0; it < maxNewton; it++) {
        assemble();
        dx.fill(0);
        cgIt += pcg(M, r, dx, nonlinear ? 1e-6 : 1e-10, 8000).it;
        iters++;
        let al = 1, rn2 = 0;
        for (let ls = 0; ls < 8; ls++) {
          for (let i = 0; i < n; i++) At[i] = dof[i] >= 0 ? A[i] + al * dx[dof[i]] : 0;
          update(At);
          rn2 = residual(r);
          if (!nonlinear || rn2 < rn * (1 - 1e-4 * al)) break;
          al *= 0.5;
        }
        A.set(At);
        rn = rn2;
        relRes = rn / fn;
        if (nonlinear && relRes < 1e-6) { converged = true; break; }
      }
    }
    // post-traitement
    const Bx = new Float64Array(nt), By = new Float64Array(nt), Bm = new Float64Array(nt), Hm = new Float64Array(nt), wE = new Float64Array(nt);
    let W = 0, lam = 0, absI = 0;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      Bx[e] = gyE[e]; By[e] = -gxE[e];
      Bm[e] = Math.hypot(Bx[e], By[e]);
      const hx = nuE[e] * (Bx[e] - p.brx), hy = nuE[e] * (By[e] - p.bry);
      Hm[e] = Math.hypot(hx, hy);
      wE[e] = p.brx || p.bry ? 0.5 * (hx * (Bx[e] - p.brx) + hy * (By[e] - p.bry)) : FT.energyDensity(p.mur, p.Js, Bm[e], nuE[e]);
      W += wE[e] * area[e];
      if (p.J) {
        const am = (A[tri[3 * e]] + A[tri[3 * e + 1]] + A[tri[3 * e + 2]]) / 3;
        lam += p.J * am * area[e];
        absI += Math.abs(p.J) * area[e];
      }
    }
    return {
      type: "mag", A, Bx, By, Bm, Hm, nu: nuE, w: wE, energy: W, linkage: lam, absI, area,
      iters, cgIt, converged, relRes, nonlinear, ms: Date.now() - T0,
    };
  };

  // props[r] = {eps (F/m), rho (C/m³), V (nombre ou null)} ; bc : 'V0' ou 'none' sur la frontière extérieure
  FT.solveElectric = function (mesh, props, bc) {
    const T0 = Date.now();
    const n = mesh.x.length, tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg;
    const { gb, gc, area } = elemGeom(mesh);
    const fixed = new Uint8Array(n), V = new Float64Array(n);
    let nfix = 0;
    if (bc === "V0") for (let i = 0; i < n; i++) if (mesh.outer[i]) { fixed[i] = 1; nfix++; }
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      if (p.V != null) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; fixed[i] = 1; V[i] = p.V; nfix++; }
    }
    let pinned = false;
    if (!nfix) { fixed[0] = 1; pinned = true; }
    const dof = new Int32Array(n);
    let nd = 0;
    for (let i = 0; i < n; i++) dof[i] = fixed[i] ? -1 : nd++;
    const M = buildStructure(tri, dof, nd);
    const f = new Float64Array(nd), va = M.va, pos = M.pos;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]], ar = area[e];
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        f[i] += (p.rho * ar) / 3;
        for (let b = 0; b < 3; b++) {
          const v = (p.eps * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
          const k = pos[9 * e + 3 * a + b];
          if (k >= 0) va[k] += v;
          else f[i] -= v * V[tri[3 * e + b]];
        }
      }
    }
    const x = new Float64Array(nd);
    const cg = pcg(M, f, x, 1e-10, 10000);
    for (let i = 0; i < n; i++) if (dof[i] >= 0) V[i] = x[dof[i]];
    const Ex = new Float64Array(nt), Ey = new Float64Array(nt), Em = new Float64Array(nt), Dm = new Float64Array(nt), wE = new Float64Array(nt);
    let W = 0;
    for (let e = 0; e < nt; e++) {
      const a0 = V[tri[3 * e]], a1 = V[tri[3 * e + 1]], a2 = V[tri[3 * e + 2]], d = 2 * area[e];
      Ex[e] = -(a0 * gb[3 * e] + a1 * gb[3 * e + 1] + a2 * gb[3 * e + 2]) / d;
      Ey[e] = -(a0 * gc[3 * e] + a1 * gc[3 * e + 1] + a2 * gc[3 * e + 2]) / d;
      Em[e] = Math.hypot(Ex[e], Ey[e]);
      const eps = props[reg[e]].eps;
      Dm[e] = eps * Em[e];
      wE[e] = 0.5 * eps * Em[e] * Em[e];
      W += wE[e] * area[e];
    }
    return { type: "elec", V, Ex, Ey, Em, Dm, w: wE, energy: W, area, cgIt: cg.it, converged: cg.res < 1e-6, pinned, ms: Date.now() - T0 };
  };

  // ---------- Régime sinusoïdal (magnétodynamique harmonique, amplitudes crête complexes) ----------
  // ∇·(ν∇A) − jωσA + σ·g_k = −Js, A = 0 sur la frontière extérieure.
  //  - bobines (Js imposé, fil fin) : courant uniforme, pas de courants induits ;
  //  - zones conductrices massives k : J = σ(−jωA + g_k), g_k (gradient de tension) inconnu,
  //    fixé par le courant total imposé I_k (0 pour une pièce non alimentée : courants induits de somme nulle).
  // Inconnues : A aux nœuds libres puis u_k = −g_k/ω (mise à l'échelle : toutes les lignes en ampères) ;
  // système complexe symétrique résolu par COCG.
  function buildStructureC(mesh, dof, nA, gOf, nG) {
    const tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg, nd = nA + nG;
    const rows = Array.from({ length: nA }, () => []);
    for (let e = 0; e < nt; e++) {
      const g = gOf[reg[e]];
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        const r = rows[i];
        for (let b = 0; b < 3; b++) { const j = dof[tri[3 * e + b]]; if (j >= 0 && r.indexOf(j) < 0) r.push(j); }
        if (g >= 0 && r.indexOf(nA + g) < 0) r.push(nA + g);
      }
    }
    const gRows = Array.from({ length: nG }, () => []);
    const stamp = new Int32Array(nA).fill(-1);
    for (let e = 0; e < nt; e++) {
      const g = gOf[reg[e]];
      if (g < 0) continue;
      for (let a = 0; a < 3; a++) { const i = dof[tri[3 * e + a]]; if (i >= 0 && stamp[i] !== g) { stamp[i] = g; gRows[g].push(i); } }
    }
    // un nœud peut être visité par la même zone après une autre : on dédoublonne
    for (let g = 0; g < nG; g++) { gRows[g] = Array.from(new Set(gRows[g])); gRows[g].push(nA + g); }
    const all = rows.concat(gRows);
    const rp = new Int32Array(nd + 1);
    for (let i = 0; i < nd; i++) rp[i + 1] = rp[i] + all[i].length;
    const ci = new Int32Array(rp[nd]), dg = new Int32Array(nd);
    for (let i = 0; i < nd; i++) {
      const r = all[i].sort((p, q) => p - q);
      for (let k = 0; k < r.length; k++) { ci[rp[i] + k] = r[k]; if (r[k] === i) dg[i] = rp[i] + k; }
    }
    const find = (i, j) => {
      let lo = rp[i], hi = rp[i + 1] - 1;
      while (lo < hi) { const m = (lo + hi) >> 1; if (ci[m] < j) lo = m + 1; else hi = m; }
      return lo;
    };
    const pos = new Int32Array(9 * nt).fill(-1), posG = new Int32Array(3 * nt).fill(-1), posGT = new Int32Array(3 * nt).fill(-1);
    for (let e = 0; e < nt; e++) {
      const g = gOf[reg[e]];
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        for (let b = 0; b < 3; b++) { const j = dof[tri[3 * e + b]]; if (j >= 0) pos[9 * e + 3 * a + b] = find(i, j); }
        if (g >= 0) { posG[3 * e + a] = find(i, nA + g); posGT[3 * e + a] = find(nA + g, i); }
      }
    }
    return { rp, ci, dg, pos, posG, posGT, vr: new Float64Array(rp[nd]), vi: new Float64Array(rp[nd]), nd };
  }

  // Gradient conjugué orthogonal complexe (COCG) préconditionné Gauss-Seidel symétrique
  function cocg(M, br, bi, xr, xi, tol, maxit) {
    const { rp, ci, vr, vi, dg, nd: n } = M;
    const rr = new Float64Array(n), ri = new Float64Array(n), zr = new Float64Array(n), zi = new Float64Array(n);
    const pr = new Float64Array(n), pi = new Float64Array(n), qr = new Float64Array(n), qi = new Float64Array(n);
    const dr = new Float64Array(n), di = new Float64Array(n);
    for (let i = 0; i < n; i++) { const a = vr[dg[i]], b = vi[dg[i]], m = a * a + b * b; dr[i] = a / m; di[i] = -b / m; }
    const mul = (ur, ui, or, oi) => {
      for (let i = 0; i < n; i++) {
        let sr = 0, si = 0;
        for (let k = rp[i]; k < rp[i + 1]; k++) { const a = vr[k], b = vi[k], c = ur[ci[k]], d = ui[ci[k]]; sr += a * c - b * d; si += a * d + b * c; }
        or[i] = sr; oi[i] = si;
      }
    };
    const prec = () => {
      for (let i = 0; i < n; i++) {
        let sr = rr[i], si = ri[i];
        for (let k = rp[i]; k < dg[i]; k++) { const a = vr[k], b = vi[k], c = zr[ci[k]], d = zi[ci[k]]; sr -= a * c - b * d; si -= a * d + b * c; }
        zr[i] = sr * dr[i] - si * di[i]; zi[i] = sr * di[i] + si * dr[i];
      }
      for (let i = n - 1; i >= 0; i--) {
        let sr = 0, si = 0;
        for (let k = dg[i] + 1; k < rp[i + 1]; k++) { const a = vr[k], b = vi[k], c = zr[ci[k]], d = zi[ci[k]]; sr += a * c - b * d; si += a * d + b * c; }
        zr[i] -= sr * dr[i] - si * di[i]; zi[i] -= sr * di[i] + si * dr[i];
      }
    };
    mul(xr, xi, qr, qi);
    let bn = 0, rn = 0;
    for (let i = 0; i < n; i++) {
      rr[i] = br[i] - qr[i]; ri[i] = bi[i] - qi[i];
      bn += br[i] * br[i] + bi[i] * bi[i]; rn += rr[i] * rr[i] + ri[i] * ri[i];
    }
    if (bn === 0) { xr.fill(0); xi.fill(0); return { it: 0, res: 0 }; }
    const target = tol * tol * bn;
    if (rn <= target) return { it: 0, res: Math.sqrt(rn / bn) };
    prec();
    pr.set(zr); pi.set(zi);
    let ror = 0, roi = 0;
    for (let i = 0; i < n; i++) { ror += rr[i] * zr[i] - ri[i] * zi[i]; roi += rr[i] * zi[i] + ri[i] * zr[i]; }
    let it = 0;
    for (; it < maxit; it++) {
      mul(pr, pi, qr, qi);
      let sr = 0, si = 0;
      for (let i = 0; i < n; i++) { sr += pr[i] * qr[i] - pi[i] * qi[i]; si += pr[i] * qi[i] + pi[i] * qr[i]; }
      const sm = sr * sr + si * si;
      if (!(sm > 0)) break;
      const ar = (ror * sr + roi * si) / sm, ai = (roi * sr - ror * si) / sm;
      rn = 0;
      for (let i = 0; i < n; i++) {
        xr[i] += ar * pr[i] - ai * pi[i]; xi[i] += ar * pi[i] + ai * pr[i];
        rr[i] -= ar * qr[i] - ai * qi[i]; ri[i] -= ar * qi[i] + ai * qr[i];
        rn += rr[i] * rr[i] + ri[i] * ri[i];
      }
      if (rn <= target) { it++; break; }
      prec();
      let r1r = 0, r1i = 0;
      for (let i = 0; i < n; i++) { r1r += rr[i] * zr[i] - ri[i] * zi[i]; r1i += rr[i] * zi[i] + ri[i] * zr[i]; }
      const om = ror * ror + roi * roi;
      const ber = (r1r * ror + r1i * roi) / om, bei = (r1i * ror - r1r * roi) / om;
      ror = r1r; roi = r1i;
      for (let i = 0; i < n; i++) {
        const tr = pr[i], ti = pi[i];
        pr[i] = zr[i] + ber * tr - bei * ti; pi[i] = zi[i] + ber * ti + bei * tr;
      }
    }
    return { it, res: Math.sqrt(rn / bn) };
  }

  // Opérations sur une matrice réelle en blocs 2×2 : structure M (rp, ci, dg), valeurs V[4k..4k+3] = [[v00, v01], [v10, v11]],
  // vecteurs entrelacés (re, im). Préconditionneur : Gauss-Seidel symétrique par blocs.
  function blockOps(M, V) {
    const { rp, ci, dg, nd: n } = M;
    const Di = new Float64Array(4 * n);
    for (let i = 0; i < n; i++) {
      const k = 4 * dg[i], a = V[k], bb = V[k + 1], c = V[k + 2], d = V[k + 3], det = a * d - bb * c;
      Di[4 * i] = d / det; Di[4 * i + 1] = -bb / det; Di[4 * i + 2] = -c / det; Di[4 * i + 3] = a / det;
    }
    const mul = (u, o) => {
      for (let i = 0; i < n; i++) {
        let s0 = 0, s1 = 0;
        for (let k = rp[i]; k < rp[i + 1]; k++) { const j = 2 * ci[k], u0 = u[j], u1 = u[j + 1]; s0 += V[4 * k] * u0 + V[4 * k + 1] * u1; s1 += V[4 * k + 2] * u0 + V[4 * k + 3] * u1; }
        o[2 * i] = s0; o[2 * i + 1] = s1;
      }
    };
    const prec = (rr, z) => {
      for (let i = 0; i < n; i++) {
        let s0 = rr[2 * i], s1 = rr[2 * i + 1];
        for (let k = rp[i]; k < dg[i]; k++) { const j = 2 * ci[k], u0 = z[j], u1 = z[j + 1]; s0 -= V[4 * k] * u0 + V[4 * k + 1] * u1; s1 -= V[4 * k + 2] * u0 + V[4 * k + 3] * u1; }
        z[2 * i] = Di[4 * i] * s0 + Di[4 * i + 1] * s1; z[2 * i + 1] = Di[4 * i + 2] * s0 + Di[4 * i + 3] * s1;
      }
      for (let i = n - 1; i >= 0; i--) {
        let s0 = 0, s1 = 0;
        for (let k = dg[i] + 1; k < rp[i + 1]; k++) { const j = 2 * ci[k], u0 = z[j], u1 = z[j + 1]; s0 += V[4 * k] * u0 + V[4 * k + 1] * u1; s1 += V[4 * k + 2] * u0 + V[4 * k + 3] * u1; }
        z[2 * i] -= Di[4 * i] * s0 + Di[4 * i + 1] * s1; z[2 * i + 1] -= Di[4 * i + 2] * s0 + Di[4 * i + 3] * s1;
      }
    };
    return { mul, prec, N: 2 * n };
  }

  // BiCGStab préconditionné : rapide mais parfois erratique ; on garde le meilleur itéré et on s'arrête s'il diverge
  function bicgstab2(M, V, b, x, tol, maxit) {
    const { mul, prec, N } = blockOps(M, V);
    const dot = (u, v) => { let q = 0; for (let i = 0; i < N; i++) q += u[i] * v[i]; return q; };
    const r = new Float64Array(N), rh = new Float64Array(N), p = new Float64Array(N), v = new Float64Array(N);
    const ph = new Float64Array(N), sv = new Float64Array(N), sh = new Float64Array(N), t = new Float64Array(N), xb = Float64Array.from(x);
    mul(x, t);
    for (let i = 0; i < N; i++) r[i] = b[i] - t[i];
    const bn = Math.sqrt(dot(b, b));
    if (bn === 0) { x.fill(0); return { it: 0, res: 0 }; }
    rh.set(r);
    let rho = 1, al = 1, om = 1, it = 0, res = Math.sqrt(dot(r, r)) / bn, best = res;
    for (; it < maxit && res > tol; it++) {
      const rho1 = dot(rh, r);
      if (!rho1) break;
      const be = (rho1 / rho) * (al / om);
      for (let i = 0; i < N; i++) p[i] = r[i] + be * (p[i] - om * v[i]);
      prec(p, ph); mul(ph, v);
      al = rho1 / dot(rh, v);
      for (let i = 0; i < N; i++) sv[i] = r[i] - al * v[i];
      const sres = Math.sqrt(dot(sv, sv)) / bn;
      if (sres <= tol) { for (let i = 0; i < N; i++) x[i] += al * ph[i]; res = sres; it++; break; }
      prec(sv, sh); mul(sh, t);
      om = dot(t, sv) / dot(t, t);
      for (let i = 0; i < N; i++) { x[i] += al * ph[i] + om * sh[i]; r[i] = sv[i] - om * t[i]; }
      rho = rho1;
      res = Math.sqrt(dot(r, r)) / bn;
      if (res < best) { best = res; xb.set(x); }
      if (!(res < 1e3 * best) || !om) break; // divergence ou blocage
    }
    if (!(res <= best)) { x.set(xb); res = best; }
    return { it, res };
  }

  // GMRES(m) préconditionné à droite : le résidu ne peut que décroître (pas de divergence)
  function gmres2(M, V, b, x, tol, maxit, m) {
    const { mul, prec, N } = blockOps(M, V);
    const dot = (u, v) => { let q = 0; for (let i = 0; i < N; i++) q += u[i] * v[i]; return q; };
    const Q = Array.from({ length: m + 1 }, () => new Float64Array(N)), Z = Array.from({ length: m }, () => new Float64Array(N));
    const H = new Float64Array((m + 1) * m), cs = new Float64Array(m), sn = new Float64Array(m), g = new Float64Array(m + 1), y = new Float64Array(m);
    const w = new Float64Array(N);
    const bn = Math.sqrt(dot(b, b));
    if (bn === 0) { x.fill(0); return { it: 0, res: 0 }; }
    let it = 0, res = 1;
    while (it < maxit) {
      mul(x, w);
      const q0 = Q[0];
      for (let i = 0; i < N; i++) q0[i] = b[i] - w[i];
      const beta = Math.sqrt(dot(q0, q0));
      res = beta / bn;
      if (res <= tol) break;
      for (let i = 0; i < N; i++) q0[i] /= beta;
      g.fill(0); g[0] = beta;
      let j = 0;
      while (j < m && it < maxit) {
        prec(Q[j], Z[j]); mul(Z[j], w);
        for (let i = 0; i <= j; i++) { const h = dot(w, Q[i]), qi = Q[i]; H[i * m + j] = h; for (let k = 0; k < N; k++) w[k] -= h * qi[k]; }
        const hn = Math.sqrt(dot(w, w));
        H[(j + 1) * m + j] = hn;
        if (hn > 0) { const qn = Q[j + 1]; for (let k = 0; k < N; k++) qn[k] = w[k] / hn; }
        for (let i = 0; i < j; i++) {
          const t = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
          H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j]; H[i * m + j] = t;
        }
        const ha = H[j * m + j], hb = H[(j + 1) * m + j], den = Math.hypot(ha, hb);
        cs[j] = den ? ha / den : 1; sn[j] = den ? hb / den : 0;
        H[j * m + j] = den; H[(j + 1) * m + j] = 0;
        g[j + 1] = -sn[j] * g[j]; g[j] = cs[j] * g[j];
        res = Math.abs(g[j + 1]) / bn;
        j++; it++;
        if (res <= tol || hn === 0) break;
      }
      for (let i = j - 1; i >= 0; i--) { let t = g[i]; for (let k = i + 1; k < j; k++) t -= H[i * m + k] * y[k]; y[i] = t / H[i * m + i]; }
      for (let i = 0; i < j; i++) { const yi = y[i], zi = Z[i]; for (let k = 0; k < N; k++) x[k] += yi * zi[k]; }
      if (res <= tol) break;
    }
    return { it, res };
  }

  // Induction crête d'un champ complexe (Bx, By) : max sur une période de |B(t)| (champ elliptique)
  const peakB = (xr, xi, yr, yi) => {
    const s = (xr * xr + xi * xi + yr * yr + yi * yi) / 2;
    const cr = xr * xr - xi * xi + yr * yr - yi * yi, cim = 2 * (xr * xi + yr * yi);
    return Math.sqrt(s + Math.hypot(cr, cim) / 2);
  };
  FT.peakB = peakB;

  // props[r] = {mur, Js, Jr, Ji (bobine, A/m²), sigma (S/m, zone massive), Ir, Ii (courant total imposé, A), solid}
  FT.solveHarmonic = function (mesh, props, omega) {
    const T0 = Date.now();
    const n = mesh.x.length, tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg, nR = mesh.nReg;
    const { gb, gc, area } = elemGeom(mesh);
    const dof = new Int32Array(n);
    let nA = 0;
    for (let i = 0; i < n; i++) dof[i] = mesh.outer[i] ? -1 : nA++;
    const gOf = new Int32Array(nR).fill(-1), gReg = [];
    for (let r = 0; r < nR; r++) if (props[r].solid) { gOf[r] = gReg.length; gReg.push(r); }
    const nG = gReg.length, nd = nA + nG;
    const M = buildStructureC(mesh, dof, nA, gOf, nG);
    const sigE = new Float64Array(nt), nuE = new Float64Array(nt), regS = new Float64Array(nR);
    let nonlinear = false;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      sigE[e] = p.solid ? p.sigma : 0;
      nuE[e] = 1 / (MU0 * p.mur);
      if (p.Js > 0) nonlinear = true;
      regS[reg[e]] += area[e];
    }
    // second membre
    const fr = new Float64Array(nd), fi = new Float64Array(nd);
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      if (!p.Jr && !p.Ji) continue;
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i >= 0) { fr[i] += (p.Jr * area[e]) / 3; fi[i] += (p.Ji * area[e]) / 3; }
      }
    }
    // ligne de contrainte : ω·C·A − jωσS·u = jI
    for (let g = 0; g < nG; g++) { const p = props[gReg[g]]; fr[nA + g] = -(p.Ii || 0); fi[nA + g] = p.Ir || 0; }
    const assemble = () => {
      const vr = M.vr, vi = M.vi, pos = M.pos;
      vr.fill(0); vi.fill(0);
      for (let e = 0; e < nt; e++) {
        const nu = nuE[e], ar = area[e], sg = sigE[e];
        for (let a = 0; a < 3; a++) {
          for (let b = 0; b < 3; b++) {
            const k = pos[9 * e + 3 * a + b];
            if (k < 0) continue;
            vr[k] += (nu * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
            if (sg) vi[k] += (omega * sg * ar * (a === b ? 2 : 1)) / 12;
          }
          if (sg && M.posG[3 * e + a] >= 0) { const c = (omega * sg * ar) / 3; vr[M.posG[3 * e + a]] += c; vr[M.posGT[3 * e + a]] += c; }
        }
      }
      for (let g = 0; g < nG; g++) { const r = gReg[g]; vi[M.dg[nA + g]] -= omega * props[r].sigma * regS[r]; }
    };
    const xr = new Float64Array(nd), xi = new Float64Array(nd);
    const Ar = new Float64Array(n), Ai = new Float64Array(n);
    const BxR = new Float64Array(nt), BxI = new Float64Array(nt), ByR = new Float64Array(nt), ByI = new Float64Array(nt), Bm = new Float64Array(nt);
    const fields = () => {
      for (let i = 0; i < n; i++) { const d = dof[i]; Ar[i] = d >= 0 ? xr[d] : 0; Ai[i] = d >= 0 ? xi[d] : 0; }
      for (let e = 0; e < nt; e++) {
        const i0 = tri[3 * e], i1 = tri[3 * e + 1], i2 = tri[3 * e + 2], d = 2 * area[e];
        const b0 = gb[3 * e], b1 = gb[3 * e + 1], b2 = gb[3 * e + 2], c0 = gc[3 * e], c1 = gc[3 * e + 1], c2 = gc[3 * e + 2];
        const gxr = (Ar[i0] * b0 + Ar[i1] * b1 + Ar[i2] * b2) / d, gxi = (Ai[i0] * b0 + Ai[i1] * b1 + Ai[i2] * b2) / d;
        const gyr = (Ar[i0] * c0 + Ar[i1] * c1 + Ar[i2] * c2) / d, gyi = (Ai[i0] * c0 + Ai[i1] * c1 + Ai[i2] * c2) / d;
        BxR[e] = gyr; BxI[e] = gyi; ByR[e] = -gxr; ByI[e] = -gxi;
        Bm[e] = peakB(gyr, gyi, -gxr, -gxi);
      }
    };
    let iters = 0, cgIt = 0, converged = true, relRes = 0;
    let any = false;
    for (let i = 0; i < nd; i++) if (fr[i] || fi[i]) { any = true; break; }
    if (any && !nonlinear) {
      assemble();
      const res = cocg(M, fr, fi, xr, xi, 1e-10, 20000);
      cgIt = res.it; iters = 1; relRes = res.res; converged = res.res < 1e-6;
      fields();
    } else if (any) {
      // Newton-Raphson sur (Re A, Im A) : ν dépend de l'induction crête, la relation n'est pas holomorphe.
      // Point fixe rapide d'abord, puis Newton (jacobien réel en blocs 2×2, GMRES) si la saturation l'exige.
      const dnuE = new Float64Array(nt), V = new Float64Array(4 * M.ci.length);
      const updNu = () => {
        for (let e = 0; e < nt; e++) {
          const p = props[reg[e]];
          if (p.Js > 0) { const q = FT.nuOfB(p.mur, p.Js, Bm[e]); nuE[e] = q[0]; dnuE[e] = q[1]; }
        }
      };
      const Rr = new Float64Array(nd), Ri = new Float64Array(nd);
      // résidu R = F − N(x) pour l'état courant (fields() et updNu() déjà appelés)
      const residual = () => {
        Rr.set(fr); Ri.set(fi);
        for (let e = 0; e < nt; e++) {
          const nu = nuE[e], ar = area[e], sg = sigE[e], g = gOf[reg[e]];
          const gxr = -ByR[e], gxi = -ByI[e], gyr = BxR[e], gyi = BxI[e];
          let sAr = 0, sAi = 0;
          for (let a = 0; a < 3; a++) { sAr += Ar[tri[3 * e + a]]; sAi += Ai[tri[3 * e + a]]; }
          for (let a = 0; a < 3; a++) {
            const i = tri[3 * e + a], d = dof[i];
            if (d < 0) continue;
            Rr[d] -= (nu * (gxr * gb[3 * e + a] + gyr * gc[3 * e + a])) / 2;
            Ri[d] -= (nu * (gxi * gb[3 * e + a] + gyi * gc[3 * e + a])) / 2;
            if (sg) {
              // jωσ ∫NaNb : Δ/12 (1+δab)
              const mr = ((sAr + Ar[i]) * ar) / 12, mi = ((sAi + Ai[i]) * ar) / 12;
              Rr[d] += omega * sg * mi; Ri[d] -= omega * sg * mr;
              if (g >= 0) { const c = (omega * sg * ar) / 3; Rr[d] -= c * xr[nA + g]; Ri[d] -= c * xi[nA + g]; }
            }
          }
          if (sg && g >= 0) { const c = (omega * sg * ar) / 3; Rr[nA + g] -= c * sAr; Ri[nA + g] -= c * sAi; }
        }
        for (let g = 0; g < nG; g++) {
          const r = gReg[g], D = omega * props[r].sigma * regS[r];
          Rr[nA + g] -= D * xi[nA + g]; Ri[nA + g] += D * xr[nA + g];
        }
        let s2 = 0;
        for (let i = 0; i < nd; i++) s2 += Rr[i] * Rr[i] + Ri[i] * Ri[i];
        return Math.sqrt(s2);
      };
      const assembleJ = (nuInit) => {
        V.fill(0);
        const pos = M.pos, qr = [0, 0, 0], qi = [0, 0, 0], wr = [0, 0, 0], wi = [0, 0, 0];
        for (let e = 0; e < nt; e++) {
          const nu = nuInit ? nuInit[e] : nuE[e], dn = nuInit ? 0 : dnuE[e], ar = area[e], sg = sigE[e];
          const gxr = -ByR[e], gxi = -ByI[e], gyr = BxR[e], gyi = BxI[e];
          for (let a = 0; a < 3; a++) { qr[a] = (gxr * gb[3 * e + a] + gyr * gc[3 * e + a]) / 2; qi[a] = (gxi * gb[3 * e + a] + gyi * gc[3 * e + a]) / 2; }
          if (dn > 0) {
            // d(Bcrête²) = Re(w·dB), w = conj(B) + ĉ·B, ĉ = conj(c)/|c|, c = Bx² + By²
            const bxr = BxR[e], bxi = BxI[e], byr = ByR[e], byi = ByI[e];
            const cr = bxr * bxr - bxi * bxi + byr * byr - byi * byi, cim = 2 * (bxr * bxi + byr * byi), cm = Math.hypot(cr, cim);
            const hr = cm > 0 ? cr / cm : 0, hi = cm > 0 ? -cim / cm : 0;
            const wxr = bxr + hr * bxr - hi * bxi, wxi = -bxi + hr * bxi + hi * bxr;
            const wyr = byr + hr * byr - hi * byi, wyi = -byi + hr * byi + hi * byr;
            for (let b = 0; b < 3; b++) { wr[b] = (wxr * gc[3 * e + b] - wyr * gb[3 * e + b]) / (2 * ar); wi[b] = (wxi * gc[3 * e + b] - wyi * gb[3 * e + b]) / (2 * ar); }
          }
          for (let a = 0; a < 3; a++) {
            for (let b = 0; b < 3; b++) {
              const k = pos[9 * e + 3 * a + b];
              if (k < 0) continue;
              const kab = (nu * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
              V[4 * k] += kab; V[4 * k + 3] += kab;
              // terme δν·qa, δν = ν'·δ(Bcrête²) = ν'·Σb Re(ωb·δAb)
              if (dn > 0) { const f = dn; V[4 * k] += f * qr[a] * wr[b]; V[4 * k + 1] -= f * qr[a] * wi[b]; V[4 * k + 2] += f * qi[a] * wr[b]; V[4 * k + 3] -= f * qi[a] * wi[b]; }
              if (sg) { const m = (omega * sg * ar * (a === b ? 2 : 1)) / 12; V[4 * k + 1] -= m; V[4 * k + 2] += m; }
            }
            if (sg && M.posG[3 * e + a] >= 0) {
              const c = (omega * sg * ar) / 3, k1 = M.posG[3 * e + a], k2 = M.posGT[3 * e + a];
              V[4 * k1] += c; V[4 * k1 + 3] += c; V[4 * k2] += c; V[4 * k2 + 3] += c;
            }
          }
        }
        for (let g = 0; g < nG; g++) {
          const r = gReg[g], D = omega * props[r].sigma * regS[r], k = M.dg[nA + g];
          V[4 * k + 1] += D; V[4 * k + 2] -= D;
        }
      };
      let fn = 0;
      for (let i = 0; i < nd; i++) fn += fr[i] * fr[i] + fi[i] * fi[i];
      fn = Math.sqrt(fn);
      converged = false;
      const b2 = new Float64Array(2 * nd), dx = new Float64Array(2 * nd), x0r = new Float64Array(nd), x0i = new Float64Array(nd);
      const setState = (ur, ui) => { xr.set(ur); xi.set(ui); fields(); updNu(); return residual(); };
      // pas de point fixe : perméabilité sécante figée, système complexe symétrique (COCG, rapide)
      const picard = () => { assemble(); cgIt += cocg(M, fr, fi, xr, xi, 1e-9, 20000).it; iters++; fields(); updNu(); return residual(); };
      // 1. point fixe tant qu'il converge vite (saturation faible ou modérée)
      let rn = fn;
      const bestR = new Float64Array(nd), bestI = new Float64Array(nd);
      let bestRn = fn, k = 0;
      for (; k < 10; k++) {
        const r2 = picard();
        if (r2 < bestRn) { bestRn = r2; bestR.set(xr); bestI.set(xi); }
        if (r2 / fn < 1e-6) { converged = true; rn = r2; break; }
        if (!(r2 < 0.3 * rn)) break; // contraction insuffisante : on passe à Newton
        rn = r2;
      }
      relRes = (converged ? rn : bestRn) / fn;
      if (!converged) {
        // 2. Newton-Raphson sur (Re A, Im A) : ν dépend de l'induction crête, relation non holomorphe.
        //    Jacobien réel en blocs 2×2 (non symétrique à cause de jωσ) : BiCGStab, GMRES en secours.
        let nuInit = null;
        if (bestRn < 0.3 * fn) rn = setState(bestR, bestI);
        else {
          // le point fixe part trop loin (forte saturation) : départ de zéro avec ν sécant à mi-saturation
          xr.fill(0); xi.fill(0); rn = setState(xr, xi);
          nuInit = Float64Array.from(nuE);
          for (let e = 0; e < nt; e++) { const p = props[reg[e]]; if (p.Js > 0) nuInit[e] = FT.nuOfB(p.mur, p.Js, 0.5 * p.Js)[0]; }
        }
        relRes = rn / fn;
        for (let it = 0; it < 60; it++) {
          assembleJ(nuInit);
          nuInit = null;
          for (let i = 0; i < nd; i++) { b2[2 * i] = Rr[i]; b2[2 * i + 1] = Ri[i]; }
          dx.fill(0);
          // Newton inexact : précision du solveur linéaire adaptée au résidu courant
          const tol = Math.min(1e-2, Math.max(1e-8, 0.1 * relRes)), rs = bicgstab2(M, V, b2, dx, tol, 3000);
          cgIt += rs.it;
          if (!(rs.res <= 10 * tol)) cgIt += gmres2(M, V, b2, dx, tol, 3000, 30).it; // repli robuste
          iters++;
          x0r.set(xr); x0i.set(xi);
          let al = 1, rn2 = rn, ok = false;
          for (let ls = 0; ls < 8; ls++) {
            for (let i = 0; i < nd; i++) { xr[i] = x0r[i] + al * dx[2 * i]; xi[i] = x0i[i] + al * dx[2 * i + 1]; }
            fields(); updNu();
            rn2 = residual();
            if (rn2 < rn * (1 - 1e-4 * al)) { ok = true; break; }
            al *= 0.5;
          }
          if (!ok) {
            // la direction de Newton n'apporte rien : un pas de point fixe, sinon on s'arrête
            setState(x0r, x0i);
            rn2 = picard();
            if (!(rn2 < rn)) { rn = setState(x0r, x0i); break; }
          }
          rn = rn2;
          relRes = rn / fn;
          if (relRes < 1e-6) { converged = true; break; }
        }
      }
    }
    fields();
    // post-traitement : densité de courant aux sommets, pertes, courant net par zone
    const gR = new Float64Array(nR), gI = new Float64Array(nR);
    for (let g = 0; g < nG; g++) { gR[gReg[g]] = -omega * xr[nA + g]; gI[gReg[g]] = -omega * xi[nA + g]; }
    const Jc = new Float64Array(3 * nt), pc = new Float64Array(3 * nt), pe = new Float64Array(nt), Hm = new Float64Array(nt), wE = new Float64Array(nt);
    const regP = new Float64Array(nR), regIr = new Float64Array(nR), regIi = new Float64Array(nR);
    const jr = [0, 0, 0], ji = [0, 0, 0];
    let W = 0, lamR = 0, lamI = 0, P = 0;
    for (let e = 0; e < nt; e++) {
      const r = reg[e], p = props[r], ar = area[e], sg = p.solid ? p.sigma : 0;
      Hm[e] = nuE[e] * Bm[e];
      wE[e] = 0.25 * nuE[e] * (BxR[e] * BxR[e] + BxI[e] * BxI[e] + ByR[e] * ByR[e] + ByI[e] * ByI[e]);
      W += wE[e] * ar;
      let mr = 0, mi = 0;
      for (let a = 0; a < 3; a++) {
        const i = tri[3 * e + a];
        mr += Ar[i] / 3; mi += Ai[i] / 3;
        if (sg) { jr[a] = sg * (omega * Ai[i] + gR[r]); ji[a] = sg * (-omega * Ar[i] + gI[r]); }
        else { jr[a] = p.Jr || 0; ji[a] = p.Ji || 0; }
        Jc[3 * e + a] = Math.hypot(jr[a], ji[a]);
        regIr[r] += (jr[a] * ar) / 3; regIi[r] += (ji[a] * ar) / 3;
      }
      if (sg) {
        // ∫|J|² sur un triangle à J linéaire = Δ/12 · Σa Σb (1+δab) Ja·Jb*
        let s2 = 0;
        for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) s2 += (a === b ? 2 : 1) * (jr[a] * jr[b] + ji[a] * ji[b]);
        const pw = (s2 * ar) / 12 / (2 * sg);
        pe[e] = pw / ar; regP[r] += pw; P += pw;
        for (let a = 0; a < 3; a++) pc[3 * e + a] = (Jc[3 * e + a] * Jc[3 * e + a]) / (2 * sg);
      }
      // flux lié des bobines : Λ = ∫ A·Js* dS
      if (!sg && (p.Jr || p.Ji)) { lamR += (mr * p.Jr + mi * p.Ji) * ar; lamI += (mi * p.Jr - mr * p.Ji) * ar; }
    }
    return {
      type: "mag", harmonic: true, omega, A: Ar, Ai, BxR, BxI, ByR, ByI, Bm, Hm, nu: nuE, w: wE, energy: W,
      Jc, pc, pe, power: P, regP, regIr, regIi, gR, gI, lamR, lamI, area,
      iters, cgIt, converged, relRes, nonlinear, ms: Date.now() - T0,
    };
  };

  // Valeur lissée aux sommets, par région (pas de lissage à travers les interfaces). Retourne 3 valeurs par triangle.
  FT.smoothPerRegion = function (mesh, val, area) {
    const nt = mesh.tri.length / 3, tri = mesh.tri, reg = mesh.reg, n = mesh.x.length;
    const out = new Float32Array(3 * nt);
    const byReg = Array.from({ length: mesh.nReg }, () => []);
    for (let e = 0; e < nt; e++) byReg[reg[e]].push(e);
    const sum = new Float64Array(n), wt = new Float64Array(n);
    for (const list of byReg) {
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; sum[i] += val[e] * area[e]; wt[i] += area[e]; }
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; out[3 * e + a] = sum[i] / wt[i]; }
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; sum[i] = 0; wt[i] = 0; }
    }
    return out;
  };

  // Lignes de niveau d'un champ nodal : retourne Float32Array [x1,y1,x2,y2,...] (mm)
  FT.contours = function (mesh, u, nLev, lo, hi) {
    let mn = Infinity, mxv = -Infinity;
    if (lo != null) { mn = lo; mxv = hi; }
    else for (let i = 0; i < u.length; i++) { if (u[i] < mn) mn = u[i]; if (u[i] > mxv) mxv = u[i]; }
    const out = [];
    if (!(mxv > mn)) return new Float32Array(0);
    const { x, y, tri } = mesh, nt = tri.length / 3, st = (mxv - mn) / nLev;
    for (let e = 0; e < nt; e++) {
      const i = tri[3 * e], j = tri[3 * e + 1], k = tri[3 * e + 2];
      const lo = Math.min(u[i], u[j], u[k]), hi = Math.max(u[i], u[j], u[k]);
      for (let l = Math.ceil((lo - mn) / st - 0.5); l <= Math.floor((hi - mn) / st - 0.5); l++) {
        const c = mn + (l + 0.5) * st;
        const pts = [];
        const cut = (a, b) => {
          if ((u[a] - c) * (u[b] - c) < 0) { const t = (c - u[a]) / (u[b] - u[a]); pts.push(x[a] + t * (x[b] - x[a]), y[a] + t * (y[b] - y[a])); }
        };
        cut(i, j); cut(j, k); cut(k, i);
        if (pts.length === 4) out.push(pts[0], pts[1], pts[2], pts[3]);
      }
    }
    return Float32Array.from(out);
  };

  // Statistiques d'une région : aire, moyenne/max du champ, énergie, courant, force (tenseur de Maxwell, méthode "coquille")
  FT.regionStats = function (mesh, res, r, depth) {
    const tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg, area = res.area;
    const fld = res.type === "mag" ? res.Bm : res.Em;
    const inBody = new Uint8Array(mesh.x.length);
    let S = 0, sum = 0, mx = 0, W = 0;
    for (let e = 0; e < nt; e++) {
      if (reg[e] !== r) continue;
      S += area[e]; sum += fld[e] * area[e]; W += res.w[e] * area[e];
      if (fld[e] > mx) mx = fld[e];
      inBody[tri[3 * e]] = inBody[tri[3 * e + 1]] = inBody[tri[3 * e + 2]] = 1;
    }
    const { gb, gc } = elemGeom(mesh);
    let fx = 0, fy = 0, contact = false;
    for (let e = 0; e < nt; e++) {
      if (reg[e] === r) continue;
      let ggx = 0, ggy = 0, any = false;
      for (let a = 0; a < 3; a++) if (inBody[tri[3 * e + a]]) { ggx += gb[3 * e + a]; ggy += gc[3 * e + a]; any = true; }
      if (!any) continue;
      let k, xx, yy, xy;
      if (res.type === "mag") {
        k = res.nu[e];
        if (res.harmonic) {
          // moyenne sur une période : <Bi·Bj> = ½ Re(Bi·Bj*)
          xx = (res.BxR[e] * res.BxR[e] + res.BxI[e] * res.BxI[e]) / 2; yy = (res.ByR[e] * res.ByR[e] + res.ByI[e] * res.ByI[e]) / 2;
          xy = (res.BxR[e] * res.ByR[e] + res.BxI[e] * res.ByI[e]) / 2;
        } else { xx = res.Bx[e] * res.Bx[e]; yy = res.By[e] * res.By[e]; xy = res.Bx[e] * res.By[e]; }
        if (Math.abs(k * FT.MU0 - 1) > 0.02) contact = true;
      } else {
        k = res.Dm[e] > 0 ? res.Dm[e] / res.Em[e] : FT.EPS0;
        xx = res.Ex[e] * res.Ex[e]; yy = res.Ey[e] * res.Ey[e]; xy = res.Ex[e] * res.Ey[e];
        if (Math.abs(k / FT.EPS0 - 1) > 0.02) contact = true;
      }
      const h2 = (xx + yy) / 2;
      // F = −∫ T·∇g dΩ, ∇g = (ggx, ggy)/(2Δ)
      fx -= (k * ((xx - h2) * ggx + xy * ggy)) / 2;
      fy -= (k * (xy * ggx + (yy - h2) * ggy)) / 2;
    }
    const o = { area: S * 1e6, mean: S ? sum / S : 0, max: mx, energy: W * depth, fx: fx * depth, fy: fy * depth, contact };
    if (res.harmonic) { o.power = res.regP[r] * depth; o.Ir = res.regIr[r]; o.Ii = res.regIi[r]; }
    return o;
  };

  // ---------- Liaison modèle -> problème ----------
  // S : {edges, labels, customMats, airMargin, meshN, ...}
  FT.buildMesh = function (S) {
    if (!S.edges.length) return null;
    const bb = FT.modelBBox(S.edges);
    const Dm = Math.max(bb.x1 - bb.x0, bb.y1 - bb.y0);
    if (!(Dm > 0)) return null;
    const m = (S.airMargin || 1) * Dm;
    const box = { x0: bb.x0 - m, y0: bb.y0 - m, x1: bb.x1 + m, y1: bb.y1 + m };
    const pslg = FT.buildPSLG(S.edges, box);
    const opt = { hmax: Dm / (S.meshN || 40), refine: 32, grad: 0.25 };
    let mesh = FT.generateMesh(pslg, opt);
    mesh.locator = FT.buildLocator(mesh);
    if (S.physics === "mag" && S.ac && S.freq > 0) {
      // sinusoïdal : mailles de δ/2 au bord des zones conductrices massives (δ = épaisseur de peau)
      FT.assignLabels(S, mesh);
      const omega = 2 * Math.PI * S.freq, dR = [];
      let dMin = Infinity;
      for (let r = 0; r < mesh.nReg; r++) {
        const l = mesh.regLabel[r], m = FT.getMaterial(S, l ? l.mat : "air"), sg = FT.sigmaOf(m);
        const stranded = l && (l.srcType === "NI" || l.srcType === "J");
        const d = sg > 0 && !stranded ? Math.sqrt(2 / (omega * MU0 * m.mur * sg)) * 1e3 : Infinity;
        dR.push(d); if (d < dMin) dMin = d;
      }
      if (dMin < 2 * opt.hmax) {
        const eps = 1e-4 * Dm, side = (x, y) => { const t = mesh.locator.find(x, y); return t < 0 ? Infinity : dR[mesh.reg[t]]; };
        opt.cap = (xa, ya, xb, yb) => {
          const L = Math.hypot(xb - xa, yb - ya) || 1, nx = (-(yb - ya) / L) * eps, ny = ((xb - xa) / L) * eps, mx = (xa + xb) / 2, my = (ya + yb) / 2;
          return 0.5 * Math.min(side(mx + nx, my + ny), side(mx - nx, my - ny));
        };
        mesh = FT.generateMesh(pslg, opt);
        mesh.locator = FT.buildLocator(mesh);
      }
      mesh.skin = { dMin, hmin: opt.hmax / opt.refine };
    }
    const g = elemGeom(mesh);
    mesh.regArea = new Float64Array(mesh.nReg); // m²
    for (let e = 0; e < mesh.reg.length; e++) mesh.regArea[mesh.reg[e]] += g.area[e];
    FT.assignLabels(S, mesh);
    return mesh;
  };

  // Associe chaque étiquette de matériau à la région qui la contient
  FT.assignLabels = function (S, mesh) {
    mesh.regLabel = new Array(mesh.nReg).fill(null);
    for (const l of S.labels) {
      const t = mesh.locator.find(l.x, l.y);
      if (t >= 0) mesh.regLabel[mesh.reg[t]] = l;
    }
  };

  FT.regionProps = function (S, mesh, physics) {
    const out = [];
    for (let r = 0; r < mesh.nReg; r++) {
      const l = mesh.regLabel[r];
      const m = FT.getMaterial(S, l ? l.mat : "air");
      if (physics === "mag") {
        let J = 0, brx = 0, bry = 0;
        const st = l ? l.srcType : "none", v = l ? +l.srcVal || 0 : 0, ph = (((l && +l.phase) || 0) * Math.PI) / 180;
        if (st === "J") J = v * 1e6;
        if (st === "NI" || st === "I") J = v / mesh.regArea[r];
        if (m.Br > 0) { const a = (((l && +l.magAngle) || 0) * Math.PI) / 180; brx = m.Br * Math.cos(a); bry = m.Br * Math.sin(a); }
        // sinusoïdal : bobine (fil fin, J uniforme) ou zone massive conductrice (courants induits, courant total imposé)
        const sigma = FT.sigmaOf(m), solid = sigma > 0 && st !== "NI" && st !== "J";
        const Iamp = st === "I" ? v : 0;
        out.push({
          mur: m.mur, Js: m.Js, J, brx, bry, sigma, solid,
          Jr: solid ? 0 : J * Math.cos(ph), Ji: solid ? 0 : J * Math.sin(ph),
          Ir: solid ? Iamp * Math.cos(ph) : 0, Ii: solid ? Iamp * Math.sin(ph) : 0,
        });
      } else {
        const fixedV = l && l.elecType === "V";
        out.push({
          eps: EPS0 * (m.cond && !fixedV ? 1e5 : m.epsr),
          rho: l && l.elecType === "rho" ? +l.elecVal || 0 : 0,
          V: fixedV ? +l.elecVal || 0 : null,
        });
      }
    }
    return out;
  };

  FT.solve = function (S, mesh) {
    FT.assignLabels(S, mesh);
    const props = FT.regionProps(S, mesh, S.physics);
    if (S.physics !== "mag") return FT.solveElectric(mesh, props, S.bcElec || "V0");
    return S.ac && S.freq > 0 ? FT.solveHarmonic(mesh, props, 2 * Math.PI * S.freq) : FT.solveMagnetic(mesh, props);
  };
})(typeof window !== "undefined" ? window : globalThis);
