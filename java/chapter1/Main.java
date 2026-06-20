package chemoinformatics;

// ============================================================
//  Main.java — Run All Chapter 1 Steps in Sequence
//  Chemoinformatics Tutorial — Chapter 1
//  Based on: Chemistry of High-Energy Materials, 4th ed.
//             Klapötke, T. M. (De Gruyter, 2019)
// ------------------------------------------------------------
//  Run in NetBeans: Right-click Main.java -> Run File (Shift+F6)
//  Or: press F6 to run the entire project
//
//  To run a SINGLE step, right-click that step's .java file
//  and choose "Run File" (Shift+F6).
// ============================================================

/**
 * Main entry point for the Chemoinformatics Tutorial — Chapter 1.
 *
 * This class runs all seven tutorial steps in sequence, printing
 * a separator between each step so you can follow along in the
 * NetBeans Output window.
 *
 * Each step corresponds exactly to a numbered step in the Python
 * notebook (Part1_Introduction_and_Setup.md):
 *
 *   Java Step 1  <->  Python: verify_install.py
 *   Java Step 2  <->  Python: step3_smiles_basics.py
 *   Java Step 3  <->  Python: step4_name_to_smiles.py
 *   Java Step 4  <->  Python: step5_canonicalize.py
 *   Java Step 5  <->  Python: step6_inchi.py
 *   Java Step 6  <->  Python: step7_draw.py
 *   Java Step 7  <->  Python: constitutional_descriptors() (Part 3 intro)
 *
 * CITATIONS:
 *   Klapötke, T.M. (2019). Chemistry of High-Energy Materials, 4th ed.
 *     De Gruyter. ISBN 978-3-11-062438-6.
 *   Willighagen, E.L. et al. (2017). CDK v2.0. J. Cheminform., 9, 33.
 *   Kim, S. et al. (2023). PubChem 2023 update. NAR 51, D1373.
 *   Weininger, D. (1988). SMILES. J. Chem. Inf. Comput. Sci., 28, 31.
 *   Heller, S.R. et al. (2015). InChI. J. Cheminform., 7, 23.
 */
public class Main {

    /** Width of the console separator line. */
    private static final int SEP_WIDTH = 62;

    public static void main(String[] args) {

        printBanner();

        // ── Step 1: Verify CDK Installation ──────────────────────────────────
        printStepHeader(1, "Verify CDK Installation");
        try {
            Step1_VerifyInstall.main(args);
        } catch (Exception e) {
            printStepError(1, e);
        }

        // ── Step 2: SMILES Basics ─────────────────────────────────────────────
        printStepHeader(2, "SMILES Basics — Parse Energetic Molecules");
        try {
            Step2_SMILESBasics.main(args);
        } catch (Exception e) {
            printStepError(2, e);
        }

        // ── Step 3: Name to SMILES (PubChem) ─────────────────────────────────
        printStepHeader(3, "Name to SMILES via PubChem REST API");
        System.out.println("  (This step requires an active internet connection.)");
        try {
            Step3_NameToSMILES.main(args);
        } catch (Exception e) {
            printStepError(3, e);
        }

        // ── Step 4: Canonical SMILES ──────────────────────────────────────────
        printStepHeader(4, "Canonical SMILES Normalization");
        try {
            Step4_Canonicalize.main(args);
        } catch (Exception e) {
            printStepError(4, e);
        }

        // ── Step 5: InChI and InChIKey ────────────────────────────────────────
        printStepHeader(5, "InChI and InChIKey Generation");
        try {
            Step5_InChI.main(args);
        } catch (Exception e) {
            printStepError(5, e);
        }

        // ── Step 6: Draw Molecules ────────────────────────────────────────────
        printStepHeader(6, "Draw Molecules to PNG Images");
        try {
            Step6_DrawMolecules.main(args);
        } catch (Exception e) {
            printStepError(6, e);
        }

        // ── Step 7: Molecular Formula and Descriptors ─────────────────────────
        printStepHeader(7, "Molecular Formula, MW, Atom Counts, Oxygen Balance");
        try {
            Step7_MolecularFormula.main(args);
        } catch (Exception e) {
            printStepError(7, e);
        }

        // ── Final summary ─────────────────────────────────────────────────────
        printFinalSummary();
    }

    // ── Helper methods ────────────────────────────────────────────────────────

    private static void printBanner() {
        String line = "=".repeat(SEP_WIDTH);
        System.out.println(line);
        System.out.println("  CHEMOINFORMATICS TUTORIAL — CHAPTER 1");
        System.out.println("  Introduction and Environment Setup");
        System.out.println("  Java / NetBeans Edition");
        System.out.println();
        System.out.println("  Reference:");
        System.out.println("    Klapötke, T.M. (2019). Chemistry of High-Energy");
        System.out.println("    Materials, 4th ed. De Gruyter.");
        System.out.println();
        System.out.println("  Library: Chemistry Development Kit (CDK) 2.9");
        System.out.println("    Willighagen et al. (2017). J. Cheminform., 9, 33.");
        System.out.println(line);
        System.out.println();
    }

    private static void printStepHeader(int stepNum, String title) {
        System.out.println();
        System.out.println("╔" + "═".repeat(SEP_WIDTH - 2) + "╗");
        System.out.printf( "║  STEP %d: %-" + (SEP_WIDTH - 12) + "s║%n", stepNum, title);
        System.out.println("╚" + "═".repeat(SEP_WIDTH - 2) + "╝");
        System.out.println();
    }

    private static void printStepError(int stepNum, Exception e) {
        System.out.println();
        System.out.printf("  [Step %d ERROR] %s: %s%n",
                          stepNum, e.getClass().getSimpleName(), e.getMessage());
        System.out.println("  Check that cdk-2.9.jar is in the lib/ folder and");
        System.out.println("  added to the NetBeans project Libraries.");
        System.out.println();
    }

    private static void printFinalSummary() {
        System.out.println();
        System.out.println("=".repeat(SEP_WIDTH));
        System.out.println("  ALL STEPS COMPLETE — CHAPTER 1 TUTORIAL");
        System.out.println("=".repeat(SEP_WIDTH));
        System.out.println();
        System.out.println("  What you have learned:");
        System.out.println("    Step 1: How to verify CDK is working in NetBeans");
        System.out.println("    Step 2: SMILES notation and molecule parsing");
        System.out.println("    Step 3: Look up compounds by name (PubChem API)");
        System.out.println("    Step 4: Canonical SMILES — one unique form per molecule");
        System.out.println("    Step 5: InChI and InChIKey — international identifiers");
        System.out.println("    Step 6: 2D structure depiction with CDK");
        System.out.println("    Step 7: Molecular formula, MW, atom counts, OB%");
        System.out.println();
        System.out.println("  Next: Chapter 2 — 3D Structures and Geometry");
        System.out.println("        (see Part2_Structures_and_3D_Geometry.md)");
        System.out.println();
        System.out.println("  Key references:");
        System.out.println("    1. Klapötke (2019). Chem. High-Energy Materials. De Gruyter.");
        System.out.println("    2. Willighagen et al. (2017). CDK v2.0. J. Cheminform. 9,33.");
        System.out.println("    3. Weininger (1988). SMILES. J. Chem. Inf. Sci. 28,31.");
        System.out.println("    4. Heller et al. (2015). InChI. J. Cheminform. 7,23.");
        System.out.println("    5. Kim et al. (2023). PubChem. NAR 51, D1373.");
        System.out.println("=".repeat(SEP_WIDTH));
    }
}
