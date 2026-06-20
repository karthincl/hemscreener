/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
/**
 *
 * @author karth
 */
package chemoinformatics.chapter2;

// ============================================================
//  Chapter 2 — Step 5: Molecular Volume, Surface Area & Density
//  Chemoinformatics & Computational Modeling of Organic Molecules
//  Reference:
//  M.Karthikeyan, Renu Vyas (2014). Practical Chemoinformatics, Springer.
//  Java 1.8 compatible — uses CDK 2.9 only
//  Run in NetBeans: Right-click -> Run File (Shift+F6)
// ============================================================
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Chapter 2, Step 5 — Molecular Volume, Surface Area and Density Estimation.
 *
 * BACKGROUND: ----------- Crystal density (ρ) is the single most important
 * physical property for predicting the performance of an explosive. Both
 * detonation velocity (D) and detonation pressure (P) scale strongly with ρ:
 *
 * D ∝ ρ (Kamlet-Jacobs, linear term) P ∝ ρ² (Kamlet-Jacobs, squared term)
 *
 * Therefore a 5% increase in density gives ~5% increase in D and ~10% increase
 * in P. Crystal density is the key design target for new energetic materials.
 * (Klapötke 2019, Chapter 1)
 *
 * METHODS FOR DENSITY PREDICTION: --------------------------------- Method 1 —
 * van der Waals volume (this step): Simple approximation using atomic vdW
 * radii. Density ≈ MW / (NA × V_vdW) Error: ±15-20% (very rough)
 *
 * Method 2 — Politzer model (recommended): ρ = α(M/V)^β + γσtot² + δν where V =
 * volume inside 0.001 au electron density isosurface (DFT), σtot² = variance of
 * ESP on the surface, ν = electrostatic balance. Error: ±2-4% (requires DFT,
 * see Chapter 6) Reference: Politzer & Murray (2016)
 *
 * Method 3 — Machine learning QSPR (Part 4 of tutorial): Train on known crystal
 * densities → predict new compounds. Error: ±3-6% with good training data.
 *
 * VAN DER WAALS RADII (Bondi 1964, used in this implementation): H: 1.20 Å C:
 * 1.70 Å N: 1.55 Å O: 1.52 Å F: 1.47 Å Cl: 1.75 Å
 *
 * GRID-BASED VOLUME ALGORITHM (this step): Place a 3D grid of points around the
 * molecule. A grid point is "inside" if it falls within the vdW radius of any
 * atom. Volume = (points inside / total points) × box volume. This is the same
 * approach used by RDKit's ComputeMolVolume().
 *
 * REFERENCES: ----------- 1. Klapötke, T.M. (2019). Chemistry of High-Energy
 * Materials, 4th ed. Chapter 1 (density-performance relationship), Ch. 2
 * (Politzer model). 2. Politzer, P.; Murray, J.S. (2016). Impact sensitivity
 * and the maximum heat of detonation. J. Mol. Model., 21, 3.
 * https://doi.org/10.1007/s00894-014-2500-3 3. Bondi, A. (1964). van der Waals
 * Volumes and Radii. J. Phys. Chem., 68(3), 441-451.
 */
public class Step5_VolumeAndDensity {

    boolean skipNullMolecules = true;
    // ── van der Waals radii in Angstrom (Bondi 1964) ─────────────────────────
    private static final Map<String, Double> VDW_RADII
            = new java.util.HashMap<String, Double>();

    static {
        VDW_RADII.put("H", 1.20);
        VDW_RADII.put("C", 1.70);
        VDW_RADII.put("N", 1.55);
        VDW_RADII.put("O", 1.52);
        VDW_RADII.put("F", 1.47);
        VDW_RADII.put("Cl", 1.75);
        VDW_RADII.put("Br", 1.85);
        VDW_RADII.put("I", 1.98);
        VDW_RADII.put("S", 1.80);
        VDW_RADII.put("P", 1.80);
        // Default for unknown elements
        VDW_RADII.put("*", 1.70);
    }

    /**
     * Return the van der Waals radius for an element symbol (Angstrom).
     */
    public static double getVdWRadius(String symbol) {
        Double r = VDW_RADII.get(symbol);
        return r != null ? r : VDW_RADII.get("*");
    }

    // ── Result container ──────────────────────────────────────────────────────
    public static class VolumeResult {

        public final String name;
        public final double avgMW;            // g/mol
        public final double vdwVolume;        // Angstrom^3
        public final double densityEstimate;  // g/cm^3 (very rough)
        public final double densityLit;       // g/cm^3 (experimental, if known)
        public final double radiusOfGyration; // Angstrom
        public final int heavyAtoms;

        public VolumeResult(String name, double mw, double vol,
                double densEst, double densLit,
                double rog, int heavyAtoms) {
            this.name = name;
            this.avgMW = mw;
            this.vdwVolume = vol;
            this.densityEstimate = densEst;
            this.densityLit = densLit;
            this.radiusOfGyration = rog;
            this.heavyAtoms = heavyAtoms;
        }
    }

    // ── Grid-based van der Waals volume ──────────────────────────────────────
    /**
     * Compute the van der Waals volume using a 3D grid method.
     *
     * Algorithm: 1. Find bounding box of all atoms + vdW margin 2. Place grid
     * points at spacing 'gridSpacing' Angstrom 3. Count grid points inside any
     * atom's vdW sphere 4. Volume = count × gridSpacing³
     *
     * This is equivalent to RDKit's AllChem.ComputeMolVolume().
     *
     * @param mol molecule with 3D coordinates
     * @param gridSpacing grid point spacing in Angstrom (0.2-0.5 typical)
     * @param boxMargin extra margin around bounding box (Angstrom)
     * @return van der Waals volume in Angstrom³
     */
    public static double computeVdWVolume(IAtomContainer mol,
            double gridSpacing,
            double boxMargin) {
        int n = mol.getAtomCount();
        if (n == 0) {
            return 0.0;
        }

        // Step 1: collect atom positions and radii
        double[] xs = new double[n], ys = new double[n],
                zs = new double[n], rs = new double[n];
        double xMin = Double.MAX_VALUE, xMax = -Double.MAX_VALUE;
        double yMin = Double.MAX_VALUE, yMax = -Double.MAX_VALUE;
        double zMin = Double.MAX_VALUE, zMax = -Double.MAX_VALUE;

        for (int i = 0; i < n; i++) {
            IAtom atom = mol.getAtom(i);
            Point3d p = atom.getPoint3d();
            if (p == null) {
                continue;
            }

            xs[i] = p.x;
            ys[i] = p.y;
            zs[i] = p.z;
            rs[i] = getVdWRadius(atom.getSymbol());

            if (p.x - rs[i] < xMin) {
                xMin = p.x - rs[i];
            }
            if (p.x + rs[i] > xMax) {
                xMax = p.x + rs[i];
            }
            if (p.y - rs[i] < yMin) {
                yMin = p.y - rs[i];
            }
            if (p.y + rs[i] > yMax) {
                yMax = p.y + rs[i];
            }
            if (p.z - rs[i] < zMin) {
                zMin = p.z - rs[i];
            }
            if (p.z + rs[i] > zMax) {
                zMax = p.z + rs[i];
            }
        }

        // Add margin
        xMin -= boxMargin;
        xMax += boxMargin;
        yMin -= boxMargin;
        yMax += boxMargin;
        zMin -= boxMargin;
        zMax += boxMargin;

        // Step 2: count interior grid points
        long inside = 0;
        long total = 0;

        for (double x = xMin; x <= xMax; x += gridSpacing) {
            for (double y = yMin; y <= yMax; y += gridSpacing) {
                for (double z = zMin; z <= zMax; z += gridSpacing) {
                    total++;
                    // Check if this grid point is inside any atom's vdW sphere
                    for (int i = 0; i < n; i++) {
                        double dx = x - xs[i];
                        double dy = y - ys[i];
                        double dz = z - zs[i];
                        if (dx * dx + dy * dy + dz * dz <= rs[i] * rs[i]) {
                            inside++;
                            break; // no double-counting
                        }
                    }
                }
            }
        }

        double cellVolume = gridSpacing * gridSpacing * gridSpacing;
        return inside * cellVolume;  // Angstrom^3
    }

    // ── Radius of gyration ────────────────────────────────────────────────────
    /**
     * Compute the radius of gyration (Rg) of a molecule. Rg measures how spread
     * out the mass is around the centroid. Compact molecules (small Rg) tend to
     * pack more efficiently in crystals.
     *
     * Rg = sqrt( Σ(mi * ri²) / Σ(mi) ) where ri is distance from atom i to the
     * centroid.
     *
     * @param mol molecule with 3D coordinates
     * @return radius of gyration in Angstrom
     */
    public static double radiusOfGyration(IAtomContainer mol) {
        int n = mol.getAtomCount();
        if (n == 0) {
            return 0.0;
        }

        // Centroid (mass-weighted; use atomic number as proxy mass)
        double totalMass = 0, cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < n; i++) {
            IAtom a = mol.getAtom(i);
            Point3d p = a.getPoint3d();
            if (p == null) {
                continue;
            }
            double mass = a.getAtomicNumber() != null ? a.getAtomicNumber() : 6.0;
            totalMass += mass;
            cx += mass * p.x;
            cy += mass * p.y;
            cz += mass * p.z;
        }
        if (totalMass == 0) {
            return 0.0;
        }
        cx /= totalMass;
        cy /= totalMass;
        cz /= totalMass;

        // Weighted mean square distance
        double sumMR2 = 0;
        for (int i = 0; i < n; i++) {
            IAtom a = mol.getAtom(i);
            Point3d p = a.getPoint3d();
            if (p == null) {
                continue;
            }
            double mass = a.getAtomicNumber() != null ? a.getAtomicNumber() : 6.0;
            double dx = p.x - cx, dy = p.y - cy, dz = p.z - cz;
            sumMR2 += mass * (dx * dx + dy * dy + dz * dz);
        }
        return Math.sqrt(sumMR2 / totalMass);
    }

    // ── Density estimation from molecular volume ──────────────────────────────
    /**
     * Rough density estimate from vdW volume.
     *
     * density (g/cm³) = MW (g/mol) / ( NA × V_mol (cm³/mol) )
     *
     * where V_mol = V_vdW × NA (converts per-molecule Å³ to per-mole cm³) and 1
     * Å³ = 1e-24 cm³
     *
     * NOTE: This systematically underestimates crystal density because: -
     * Packing efficiency in crystals is ~65-75% (not 100%) - The method ignores
     * intermolecular interactions The Politzer DFT method (Chapter 6) gives
     * much better accuracy.
     *
     * @param mw_g_mol molecular weight (g/mol)
     * @param vdw_ang3 van der Waals volume (Angstrom^3 per molecule)
     * @return estimated density in g/cm^3
     */
    public static double estimateDensity(double mw_g_mol, double vdw_ang3) {
        final double NA = 6.02214076e23;  // Avogadro's number
        final double ANG3_TO_CM3 = 1.0e-24;       // 1 Å³ = 1e-24 cm³
        double V_molar_cm3 = vdw_ang3 * ANG3_TO_CM3 * NA;
        return mw_g_mol / V_molar_cm3;
    }

    // ── Full analysis for one compound ────────────────────────────────────────
    public static VolumeResult analyseCompound(
            String name, String smiles, double litDensity,
            SmilesParser parser) {

        // Generate 3D structure
        Step1_3DGeneration.Structure3D s
                = Step1_3DGeneration.generateStructure3D(name, smiles, parser);

        if (!s.success || s.mol == null) {
            System.out.println("  3D generation failed for " + name + ": " + s.message);
            return new VolumeResult(name, 0, 0, 0, litDensity, 0, 0);
        }

        IAtomContainer mol = s.mol;

        // Molecular weight
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        IMolecularFormula mf
                = MolecularFormulaManipulator.getMolecularFormula(mol);
        double mw = MolecularFormulaManipulator.getTotalNaturalAbundance(mf);

        // van der Waals volume (grid spacing 0.3 Å — balance of speed and accuracy)
        double vdwVol = computeVdWVolume(mol, 0.3, 2.0);

        // Density estimate
        double densEst = estimateDensity(mw, vdwVol);

        // Radius of gyration
        double rog = radiusOfGyration(mol);

        // Heavy atom count
        int heavyAtoms = 0;
        for (int i = 0; i < mol.getAtomCount(); i++) {
            if (!mol.getAtom(i).getSymbol().equals("H")) {
                heavyAtoms++;
            }
        }

        return new VolumeResult(name, mw, vdwVol, densEst,
                litDensity, rog, heavyAtoms);
    }

    // ── Main ─────────────────────────────────────────────────────────────────
    public static void main(String[] args) {

        System.out.println(rep('=', 70));
        System.out.println("  CHAPTER 2 - STEP 5: Molecular Volume, Surface Area & Density");
        System.out.println("  van der Waals Volume | Radius of Gyration | Density Estimate");
        System.out.println("  Based on: Klapötke (2019); Politzer & Murray (2016)");
        System.out.println(rep('=', 70));
        System.out.println();

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // ── Dataset with literature crystal densities ─────────────────────────
        // Experimental crystal densities from Klapötke (2019), Table 1.2
        // and Dobratz (1981), LLNL Explosives Handbook
        Object[][] dataset = {
            // {name,  SMILES,  lit_density_g_cm3}
            {"Nitromethane", "C[N+](=O)[O-]", 1.137},
            {"TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", 1.654},
            {"Picric acid", "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", 1.763},
            {"RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]", 1.806},
            {"HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1", 1.906},
            {"PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]", 1.778},
            {"FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N", 1.878},
            {"TATB", "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]", 1.938},};

        System.out.printf("  %-18s  %8s  %10s  %10s  %10s  %8s%n",
                "Compound", "MW", "vdW Vol", "Est.Dens.", "Lit.Dens.", "Rg (Å)");
        System.out.printf("  %-18s  %8s  %10s  %10s  %10s  %8s%n",
                "", "(g/mol)", "(Å³)", "(g/cm³)", "(g/cm³)", "");
        System.out.println("  " + rep('-', 75));

        List<VolumeResult> results = new ArrayList<VolumeResult>();

        for (Object[] row : dataset) {
            String name = (String) row[0];
            String smiles = (String) row[1];
            double litDens = (Double) row[2];

            VolumeResult r = analyseCompound(name, smiles, litDens, parser);
            results.add(r);

            double pctErr = (r.densityLit > 0 && r.densityEstimate > 0)
                    ? 100.0 * (r.densityEstimate - r.densityLit) / r.densityLit
                    : 0;

            System.out.printf("  %-18s  %8.3f  %10.2f  %10.3f  %10.3f  %8.3f%n",
                    r.name, r.avgMW, r.vdwVolume,
                    r.densityEstimate, r.densityLit, r.radiusOfGyration);
        }

        // ── Density comparison summary ─────────────────────────────────────────
        System.out.println();
        System.out.println("  DENSITY ESTIMATION ERROR ANALYSIS:");
        System.out.println("  " + rep('-', 55));
        System.out.printf("  %-18s  %12s  %10s  %8s%n",
                "Compound", "Est (g/cm³)", "Lit (g/cm³)", "Error%");
        System.out.println("  " + rep('-', 55));

        double sumAbsErr = 0;
        int countErr = 0;

        for (VolumeResult r : results) {
            if (r.densityEstimate > 0 && r.densityLit > 0) {
                double err = 100.0 * (r.densityEstimate - r.densityLit) / r.densityLit;
                sumAbsErr += Math.abs(err);
                countErr++;
                System.out.printf("  %-18s  %12.3f  %10.3f  %+8.1f%%%n",
                        r.name, r.densityEstimate, r.densityLit, err);
            }
        }
        if (countErr > 0) {
            System.out.printf("  %-18s  %12s  %10s  %8.1f%%%n",
                    "Mean abs error", "", "",
                    sumAbsErr / countErr);
        }

        // ── Method comparison and recommendations ─────────────────────────────
        System.out.println();
        System.out.println("  DENSITY PREDICTION METHOD COMPARISON:");
        System.out.println("  " + rep('-', 62));
        System.out.printf("  %-28s  %10s  %s%n", "Method", "Error (%)", "Requires");
        System.out.println("  " + rep('-', 62));
        System.out.printf("  %-28s  %10s  %s%n",
                "vdW volume (this step)", "15-20%", "3D coords only");
        System.out.printf("  %-28s  %10s  %s%n",
                "QSPR (ML model, Part 4)", "3-6%", "descriptor set");
        System.out.printf("  %-28s  %10s  %s%n",
                "Politzer DFT (Chapter 6)", "2-4%", "B3LYP/6-31G(d)");
        System.out.printf("  %-28s  %10s  %s%n",
                "X-ray crystallography", "<1%", "synthesized crystal");
        System.out.println();
        System.out.println("  CONCLUSION:");
        System.out.println("  The vdW volume method is useful for rapid structural comparison");
        System.out.println("  and ranking, but not for accurate density prediction.");
        System.out.println("  For publication-quality density values, use the Politzer model");
        System.out.println("  (requires DFT calculation — see Chapter 6).");
        System.out.println("  Reference: Politzer & Murray (2016). J. Mol. Model. 21, 3.");
        System.out.println(rep('=', 70));
    }

    static String rep(char c, int n) {
        return new String(new char[n]).replace('\0', c);
    }
}
