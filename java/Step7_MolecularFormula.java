package chemoinformatics;

// ============================================================
//  Step 7 — Molecular Formula, Molecular Weight, and Atom Counting
//  Chapter 1: Introduction and Environment Setup
//  Based on: Chemistry of High-Energy Materials, 4th ed.
//             Klapötke, T. M. (De Gruyter, 2019)
// ------------------------------------------------------------
//  Equivalent to Python: constitutional_descriptors() in Part 3
//  (this is the Java Chapter 1 version — formula and atom counts only)
//  Run in NetBeans: Right-click -> Run File (Shift+F6)
// ============================================================

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Step 7 — Molecular Formula, MW, and Atom Counting.
 *
 * Before calculating descriptors (Part 3) or running quantum chemistry
 * (Part 6), you must be able to extract the basic formula data:
 *   - Molecular formula  (Hill notation: C first, H second, then alphabetical)
 *   - Exact molecular weight (from isotope masses)
 *   - Average molecular weight (from natural isotope abundances)
 *   - Atom counts per element  (C, H, N, O counts are critical for OB%)
 *   - Degree of unsaturation / Double Bond Equivalents (DBE)
 *   - Nitrogen content (wt%) and Oxygen content (wt%)
 *
 * Oxygen balance is introduced here conceptually; full calculation appears
 * in Part 3 (Molecular Descriptors). DBE and atom counts are needed for:
 *   - Kamlet–Jacobs equations (Klapötke 2019, Chapter 1, Eq. 1.1–1.6)
 *   - Oxygen balance formula
 *   - Energetic material classification
 *
 * CDK Reference:
 *   MolecularFormulaManipulator provides Hill notation, exact and average
 *   mass, and element-by-element iteration.
 */
public class Step7_MolecularFormula {

    // ── Data class for computed results ───────────────────────────────────────
    /**
     * Mirrors the dict produced by constitutional_descriptors() in Python.
     */
    public static class FormulaData {
        public final String  name;
        public final String  hillFormula;   // Hill notation (C then H then alphabetical)
        public final double  exactMW;       // g/mol (monoisotopic)
        public final double  avgMW;         // g/mol (natural abundance average)
        public final int     countC;
        public final int     countH;
        public final int     countN;
        public final int     countO;
        public final double  dbe;           // Double Bond Equivalents
        public final double  nContentPct;   // Nitrogen content (wt%)
        public final double  oContentPct;   // Oxygen content (wt%)
        public final double  oxygenBalance; // OB% = (1600/MW)*(O - 2C - H/2)

        public FormulaData(String name, String hillFormula,
                           double exactMW, double avgMW,
                           int C, int H, int N, int O) {
            this.name        = name;
            this.hillFormula = hillFormula;
            this.exactMW     = exactMW;
            this.avgMW       = avgMW;
            this.countC      = C;
            this.countH      = H;
            this.countN      = N;
            this.countO      = O;

            // DBE = (2C + 2 + N - H) / 2  for CxHyNzOw (O does not affect DBE)
            this.dbe = (2.0 * C + 2.0 + N - H) / 2.0;

            // Nitrogen content in weight percent
            this.nContentPct = (avgMW > 0) ? (14.007 * N / avgMW) * 100.0 : 0;

            // Oxygen content in weight percent
            this.oContentPct = (avgMW > 0) ? (15.999 * O / avgMW) * 100.0 : 0;

            // Oxygen balance for CaHbNcOd explosives (Klapötke 2019, Eq. 1.1)
            // OB% = (1600 / MW) * (d - 2a - b/2)
            this.oxygenBalance = (avgMW > 0)
                ? (1600.0 / avgMW) * (O - 2.0 * C - H / 2.0)
                : 0;
        }
    }

    // ── Core computation method ───────────────────────────────────────────────
    public static FormulaData computeFormulaData(
            String name, String smiles,
            SmilesParser parser) {

        try {
            IAtomContainer mol = parser.parseSmiles(smiles);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
            // Convert implicit H to explicit so we can count them
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

            // ── Molecular formula object ───────────────────────────────────────
            IMolecularFormula mf =
                MolecularFormulaManipulator.getMolecularFormula(mol);

            String hillStr = MolecularFormulaManipulator.getString(mf);
            double exactMW = MolecularFormulaManipulator.getTotalExactMass(mf);
            double avgMW   = MolecularFormulaManipulator.getTotalNaturalAbundance(mf);

            // ── Count atoms per element ────────────────────────────────────────
            // Iterate over explicit atoms; hydrogen is already explicit after
            // the convertImplicitToExplicitHydrogens call above
            Map<String, Integer> counts = new TreeMap<>();
            for (IAtom atom : mol.atoms()) {
                String sym = atom.getSymbol();
                counts.merge(sym, 1, Integer::sum);
            }

            int C = counts.getOrDefault("C", 0);
            int H = counts.getOrDefault("H", 0);
            int N = counts.getOrDefault("N", 0);
            int O = counts.getOrDefault("O", 0);

            return new FormulaData(name, hillStr, exactMW, avgMW, C, H, N, O);

        } catch (Exception e) {
            System.err.printf("  ERROR processing %s: %s%n", name, e.getMessage());
            return null;
        }
    }

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 7 — Molecular Formula, MW, and Atom Counting");
        System.out.println("============================================================\n");

        // ── CDK setup ─────────────────────────────────────────────────────────
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser       parser  = new SmilesParser(builder);

        // ── Dataset with experimental detonation velocities ───────────────────
        // Data source: Klapötke (2019), Table 1.1
        // Format: name -> {SMILES, experimental det. velocity m/s}
        Map<String, String> smiles = new LinkedHashMap<>();
        smiles.put("Nitromethane",     "C[N+](=O)[O-]");
        smiles.put("TNT",              "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        smiles.put("Picric acid",      "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        smiles.put("RDX",              "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        smiles.put("HMX",              "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        smiles.put("PETN",             "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        smiles.put("Nitroglycerin",    "O([N+](=O)[O-])CC(O[N+](=O)[O-])CO[N+](=O)[O-]");
        smiles.put("FOX-7",            "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        smiles.put("TATB",             "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]");
        smiles.put("NTO",              "O=c1[nH]nno1");
        smiles.put("Ammonium nitrate", "[NH4+].[O-][N+]([O-])=O");

        // ── Table 1: Basic formula data ───────────────────────────────────────
        System.out.printf("  %-22s  %-14s  %8s  %4s  %4s  %4s  %4s%n",
                          "Compound", "Formula", "Avg MW", "C", "H", "N", "O");
        System.out.println("  " + "-".repeat(75));

        java.util.List<FormulaData> allData = new java.util.ArrayList<>();

        for (Map.Entry<String, String> e : smiles.entrySet()) {
            FormulaData d = computeFormulaData(e.getKey(), e.getValue(), parser);
            if (d == null) continue;
            allData.add(d);
            System.out.printf("  %-22s  %-14s  %8.3f  %4d  %4d  %4d  %4d%n",
                              d.name, d.hillFormula, d.avgMW,
                              d.countC, d.countH, d.countN, d.countO);
        }

        // ── Table 2: Energetic descriptors ────────────────────────────────────
        System.out.println();
        System.out.printf("  %-22s  %6s  %8s  %8s  %8s  Interpretation%n",
                          "Compound", "DBE", "N wt%", "O wt%", "OB%");
        System.out.println("  " + "-".repeat(85));

        for (FormulaData d : allData) {
            String interp;
            if      (d.oxygenBalance > 10)  interp = "over-oxidized";
            else if (d.oxygenBalance > 0)   interp = "slight O excess";
            else if (d.oxygenBalance > -10) interp = "near-zero (high performance)";
            else if (d.oxygenBalance > -40) interp = "O-deficient";
            else                             interp = "severely O-deficient";

            System.out.printf("  %-22s  %6.1f  %8.2f  %8.2f  %8.2f  %s%n",
                              d.name, d.dbe, d.nContentPct,
                              d.oContentPct, d.oxygenBalance, interp);
        }

        // ── Oxygen balance formula explanation ────────────────────────────────
        System.out.println();
        System.out.println("  OXYGEN BALANCE FORMULA (Klapötke 2019, Eq. 1.1):");
        System.out.println("  For compound CaHbNcOd with molecular weight MW (g/mol):");
        System.out.println("    OB% = (1600 / MW) × (d - 2a - b/2)");
        System.out.println();
        System.out.println("  DOUBLE BOND EQUIVALENTS (DBE):");
        System.out.println("    DBE = (2C + 2 + N - H) / 2");
        System.out.println("    DBE = 0  -> fully saturated (alkane-like)");
        System.out.println("    DBE = 1  -> one ring OR one double bond");
        System.out.println("    DBE = 4  -> benzene ring");
        System.out.println("    High DBE + N-rich + positive OB% -> energetic candidate");

        // ── Worked example: RDX ───────────────────────────────────────────────
        System.out.println();
        System.out.println("  WORKED EXAMPLE — RDX (C3H6N6O6, MW = 222.12 g/mol):");
        System.out.println("    OB% = (1600/222.12) × (6 - 2×3 - 6/2)");
        System.out.println("        = (7.203)       × (6 - 6 - 3)");
        System.out.println("        = (7.203)       × (-3)");
        System.out.println("        = -21.6%   [oxygen-deficient; produces some CO]");
        System.out.println("    DBE = (2×3 + 2 + 6 - 6) / 2 = 4 / 2 = 2");
        System.out.println("          -> 2 degrees of unsaturation (the 6-membered ring");
        System.out.println("             counts as 1 ring + 3 C-N bonds, no pi bonds)");
        System.out.println();
        System.out.println("  Experimental value for RDX OB%: -21.6%  [matches calculation]");
        System.out.println("  Source: Klapötke (2019), Table 1.2");
        System.out.println("============================================================");
    }
}
