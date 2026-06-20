package chemoinformatics;

// ============================================================
//  Step 6 — Draw Molecules (2D Structure Images)
//  Chapter 1: Introduction and Environment Setup
//  Based on: Chemistry of High-Energy Materials, 4th ed.
//             Klapötke, T. M. (De Gruyter, 2019)
// ------------------------------------------------------------
//  Equivalent to Python:  step7_draw.py  (uses RDKit Draw)
//  Java equivalent uses CDK's DepictionGenerator (CDK 2.x)
//
//  Output: PNG image files saved to molecules/ subfolder
//  Run in NetBeans: Right-click -> Run File (Shift+F6)
// ============================================================

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.renderer.color.UniColor;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Step 6 — Draw Molecules to PNG Image Files.
 *
 * CDK's DepictionGenerator produces publication-quality 2D structure
 * diagrams. The StructureDiagramGenerator computes 2D coordinates from
 * a SMILES string (which has no geometry information).
 *
 * Output directory: molecules/   (created automatically in project folder)
 *   individual_structures/TNT.png, RDX.png, ...  — one PNG per compound
 *   grid/energetic_molecules_grid.png             — 3x3 grid image
 *
 * Python equivalent:
 *   Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(400,300))
 *
 * CDK DepictionGenerator reference:
 *   Hanson, R.M. et al. (2017). The CDK depiction module.
 *   Part of the CDK 2.0 release. https://cdk.github.io
 *
 * Energetic compound structures from:
 *   Klapötke, T.M. (2019). Chemistry of High-Energy Materials, 4th ed.
 *   Figures 1.1–1.8. De Gruyter.
 */
public class Step6_DrawMolecules {

    // ── Energetic compound dataset ────────────────────────────────────────────
    // Listed in order of increasing complexity (Klapötke 2019, Chapter 1)
    private static final Map<String, String> COMPOUNDS = new LinkedHashMap<>();
    static {
        COMPOUNDS.put("Nitromethane",  "C[N+](=O)[O-]");
        COMPOUNDS.put("TNT",           "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        COMPOUNDS.put("Picric acid",   "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        COMPOUNDS.put("RDX",           "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        COMPOUNDS.put("HMX",           "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        COMPOUNDS.put("PETN",          "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        COMPOUNDS.put("FOX-7",         "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        COMPOUNDS.put("TATB",          "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]");
        COMPOUNDS.put("NTO",           "O=c1[nH]nno1");
    }

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 6 — Draw Molecules to PNG Images");
        System.out.println("  CDK DepictionGenerator + StructureDiagramGenerator");
        System.out.println("============================================================\n");

        // ── Create output directories ─────────────────────────────────────────
        File outDir     = new File("molecules");
        File indivDir   = new File("molecules/individual_structures");
        File gridDir    = new File("molecules/grid");

        outDir.mkdirs();
        indivDir.mkdirs();
        gridDir.mkdirs();

        // ── Setup CDK parser and diagram generator ────────────────────────────
        IChemObjectBuilder    builder = SilentChemObjectBuilder.getInstance();
        SmilesParser          parser  = new SmilesParser(builder);
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();

        // Configure the depiction generator
        // White background, standard bond width, atom coloring by CPK convention
        DepictionGenerator dg = new DepictionGenerator()
            .withAtomColors()           // CPK colors: C=black, N=blue, O=red, H=grey
            .withBackgroundColor(Color.WHITE)
            .withPadding(0.5)           // padding around structure (bond-length units)
            .withTerminalCarbons()      // show terminal CH3 carbon labels
            .withAromaticDisplay();     // draw aromatic rings with circle notation

        // ── Process each compound ─────────────────────────────────────────────
        List<IAtomContainer> molList    = new ArrayList<>();
        List<String>         labelList  = new ArrayList<>();
        int successCount = 0;

        System.out.printf("  %-22s  %-10s  %s%n", "Compound", "Status", "Output file");
        System.out.println("  " + "-".repeat(65));

        for (Map.Entry<String, String> entry : COMPOUNDS.entrySet()) {
            String name   = entry.getKey();
            String smiles = entry.getValue();
            String safeN  = name.replace(" ", "_").replace("-", "_");

            try {
                // 1. Parse SMILES -> molecule graph
                IAtomContainer mol = parser.parseSmiles(smiles);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);

                // 2. Generate 2D coordinates
                //    SMILES has no geometry; SDG computes a layout using
                //    template-matching and ring-system rules
                sdg.setMolecule(mol);
                sdg.generateCoordinates();
                IAtomContainer mol2D = sdg.getMolecule();

                // 3. Draw individual PNG (400 x 300 pixels)
                File outFile = new File(indivDir, safeN + ".png");
                dg.withSize(400, 300)
                  .depict(mol2D)
                  .writeTo(outFile.getAbsolutePath());

                System.out.printf("  %-22s  %-10s  %s%n",
                                  name, "OK", outFile.getPath());

                // Collect for grid image
                molList.add(mol2D);
                labelList.add(name);
                successCount++;

            } catch (Exception e) {
                System.out.printf("  %-22s  %-10s  ERROR: %s%n",
                                  name, "FAILED", e.getMessage());
            }
        }

        // ── Draw grid image ───────────────────────────────────────────────────
        System.out.println();
        if (molList.size() >= 2) {
            try {
                File gridFile = new File(gridDir, "energetic_molecules_grid.png");

                // Convert list to array for CDK API
                IAtomContainer[] molArray    = molList.toArray(new IAtomContainer[0]);
                String[]         labelArray  = labelList.toArray(new String[0]);

                // CDK grid depiction: 3 columns, auto rows
                dg.withSize(1200, 900)             // total image size
                  .depict(molList, 3, -1)          // 3 columns, auto rows
                  .writeTo(gridFile.getAbsolutePath());

                System.out.printf("  Grid image saved : %s%n", gridFile.getPath());

            } catch (Exception e) {
                System.out.println("  WARNING: Grid image generation failed: " + e.getMessage());
                System.out.println("  Individual PNG files are still available.");
            }
        }

        // ── Display instructions for NetBeans ─────────────────────────────────
        System.out.println();
        System.out.println("  HOW TO VIEW IN NETBEANS:");
        System.out.println("  1. In the Projects panel, right-click the project root");
        System.out.println("  2. Select 'Refresh' to see the new 'molecules' folder");
        System.out.println("  3. Double-click any .png file to open it in the viewer");
        System.out.println("  4. Or: File menu -> Open File -> navigate to molecules/");
        System.out.println();
        System.out.printf("  Processed: %d/%d compounds successfully%n",
                          successCount, COMPOUNDS.size());
        System.out.println();
        System.out.println("  Color scheme (CPK convention):");
        System.out.println("    Carbon   (C) : black");
        System.out.println("    Nitrogen (N) : dark blue");
        System.out.println("    Oxygen   (O) : red");
        System.out.println("    Hydrogen (H) : light grey (shown only on heteroatoms)");
        System.out.println();
        System.out.println("  Depiction reference:");
        System.out.println("    CDK DepictionGenerator — https://cdk.github.io");
        System.out.println("    Structures from Klapötke (2019), Figs. 1.1-1.8");
        System.out.println("============================================================");
    }
}
