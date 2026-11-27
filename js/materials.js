// Bibliothèque de matériaux + modèle de saturation B(H)
(function (G) {
  "use strict";
  const FT = (G.FT = G.FT || {});
  const MU0 = 4e-7 * Math.PI;
  FT.MU0 = MU0;
  FT.EPS0 = 8.8541878128e-12;

  // mur : perméabilité relative initiale ; Js : polarisation à saturation (T), 0 = linéaire
  // Br : rémanence (T) pour les aimants ; epsr : permittivité relative ; cond : conducteur (électrostatique)
  // sigma : conductivité (MS/m) pour le régime sinusoïdal ; 0 pour les tôles feuilletées et les isolants
  FT.defaultMaterials = [
    { id: "air", name: "Air", color: "#ffffff", mur: 1, Js: 0, Br: 0, epsr: 1, cond: false, sigma: 0 },
    { id: "cu", name: "Cuivre", color: "#e9a45b", mur: 1, Js: 0, Br: 0, epsr: 1, cond: true, sigma: 58 },
    { id: "alu", name: "Aluminium", color: "#c9d3dc", mur: 1, Js: 0, Br: 0, epsr: 1, cond: true, sigma: 35 },
    { id: "fe_lin", name: "Fer doux linéaire (µr 2000)", color: "#8fa3b5", mur: 2000, Js: 0, Br: 0, epsr: 1, cond: true, sigma: 10 },
    { id: "steel", name: "Acier doux (saturable, Js 2,05 T)", color: "#7f93a8", mur: 1500, Js: 2.05, Br: 0, epsr: 1, cond: true, sigma: 6 },
    { id: "fesi", name: "Tôle FeSi feuilletée (saturable, Js 1,95 T)", color: "#9aa9c9", mur: 5000, Js: 1.95, Br: 0, epsr: 1, cond: true, sigma: 0 },
    { id: "ferrite", name: "Ferrite MnZn (saturable, Js 0,45 T)", color: "#77808a", mur: 2200, Js: 0.45, Br: 0, epsr: 1, cond: false, sigma: 0 },
    { id: "powder", name: "Poudre de fer (µr 60, Js 1,2 T)", color: "#b3a58c", mur: 60, Js: 1.2, Br: 0, epsr: 1, cond: false, sigma: 0 },
    { id: "ndfeb", name: "Aimant NdFeB N42 (Br 1,30 T)", color: "#e27a7a", mur: 1.05, Js: 0, Br: 1.3, epsr: 1, cond: true, sigma: 0.67 },
    { id: "smco", name: "Aimant SmCo (Br 1,05 T)", color: "#d98fb8", mur: 1.05, Js: 0, Br: 1.05, epsr: 1, cond: true, sigma: 1.2 },
    { id: "ferr_mag", name: "Aimant ferrite (Br 0,40 T)", color: "#b98a8a", mur: 1.1, Js: 0, Br: 0.4, epsr: 1, cond: false, sigma: 0 },
    { id: "fr4", name: "FR4 (εr 4,4)", color: "#b9d8a6", mur: 1, Js: 0, Br: 0, epsr: 4.4, cond: false, sigma: 0 },
    { id: "ptfe", name: "PTFE (εr 2,1)", color: "#e4e9c8", mur: 1, Js: 0, Br: 0, epsr: 2.1, cond: false, sigma: 0 },
    { id: "kapton", name: "Kapton (εr 3,4)", color: "#e8c877", mur: 1, Js: 0, Br: 0, epsr: 3.4, cond: false, sigma: 0 },
    { id: "alumina", name: "Alumine (εr 9,8)", color: "#d9d2e8", mur: 1, Js: 0, Br: 0, epsr: 9.8, cond: false, sigma: 0 },
    { id: "water", name: "Eau (εr 80)", color: "#a9d4ef", mur: 1, Js: 0, Br: 0, epsr: 80, cond: false, sigma: 0 },
  ];

  FT.sigmaOf = (m) => Math.max(0, +(m && m.sigma) || 0) * 1e6; // S/m

  FT.getMaterial = function (S, id) {
    const all = FT.defaultMaterials.concat((S && S.customMats) || []);
    for (const m of all) if (m.id === id) return m;
    return FT.defaultMaterials[0];
  };

  // Modèle : B = µ0·H + Js·(2/π)·atan(π(µr−1)µ0·H / (2Js))
  // Retourne [ν, dν/d(B²)] avec ν = H/B
  FT.nuOfB = function (mur, Js, B) {
    const nu0 = 1 / (MU0 * mur);
    if (!(Js > 0) || B < 1e-7) return [nu0, 0];
    const k = (Math.PI * (mur - 1) * MU0) / (2 * Js), c = (Js * 2) / Math.PI;
    let H = B * nu0;
    for (let it = 0; it < 80; it++) {
      const f = MU0 * H + c * Math.atan(k * H) - B;
      const d = MU0 + (c * k) / (1 + k * k * H * H);
      const dH = -f / d;
      H += dH;
      if (Math.abs(dH) < 1e-11 * Math.abs(H) + 1e-12) break;
    }
    const d = MU0 + (c * k) / (1 + k * k * H * H);
    const nu = H / B;
    return [nu, Math.max(0, (1 / d - nu) / (2 * B * B))];
  };

  // Densité d'énergie magnétique ∫H·dB (J/m³)
  FT.energyDensity = function (mur, Js, B, nu) {
    if (!(Js > 0)) return 0.5 * nu * B * B;
    const H = nu * B;
    const k = (Math.PI * (mur - 1) * MU0) / (2 * Js), c = (Js * 2) / Math.PI;
    const co = (MU0 * H * H) / 2 + c * (H * Math.atan(k * H) - Math.log(1 + k * k * H * H) / (2 * k));
    return B * H - co;
  };
})(typeof window !== "undefined" ? window : globalThis);
