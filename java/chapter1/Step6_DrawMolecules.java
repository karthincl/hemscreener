package chemoinformatics;

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

public class Step6_DrawMolecules {

    // ── Energetic compound dataset ────────────────────────────────────────────
    // Listed in order of increasing complexity (Klapötke 2019, Chapter 1)
    private static final Map<String, String> COMPOUNDS = new LinkedHashMap<>();

    static {
        COMPOUNDS.put("Nitromethane", "C[N+](=O)[O-]");
        COMPOUNDS.put("TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        COMPOUNDS.put("Picric acid", "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        COMPOUNDS.put("RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        COMPOUNDS.put("HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        COMPOUNDS.put("PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        COMPOUNDS.put("FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        COMPOUNDS.put("TATB", "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]");
        COMPOUNDS.put("NTO", "O=c1[nH]nno1");
    }

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 6 — Draw Molecules to PNG Images");
        System.out.println("  CDK DepictionGenerator + StructureDiagramGenerator");
        System.out.println("============================================================\n");

        // ── Create output directories ─────────────────────────────────────────
        File outDir = new File("molecules");
        File indivDir = new File("molecules/individual_structures");
        File gridDir = new File("molecules/grid");

        outDir.mkdirs();
        indivDir.mkdirs();
        gridDir.mkdirs();

        // ── Setup CDK parser and diagram generator ────────────────────────────
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();

        DepictionGenerator dg = new DepictionGenerator()
                .withAtomColors() // CPK colors: C=black, N=blue, O=red, H=grey
                .withBackgroundColor(Color.WHITE)
                .withPadding(0.5) // padding around structure (bond-length units)
                .withTerminalCarbons() // show terminal CH3 carbon labels
                .withAromaticDisplay();     // draw aromatic rings with circle notation

        // ── Process each compound ─────────────────────────────────────────────
        List<IAtomContainer> molList = new ArrayList<>();
        List<String> labelList = new ArrayList<>();
        int successCount = 0;

        System.out.printf("  %-22s  %-10s  %s%n", "Compound", "Status", "Output file");
        System.out.println("  " + "-".repeat(65));

        for (Map.Entry<String, String> entry : COMPOUNDS.entrySet()) {
            String name = entry.getKey();
            String smiles = entry.getValue();
            String safeN = name.replace(" ", "_").replace("-", "_");

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
                IAtomContainer[] molArray = molList.toArray(new IAtomContainer[0]);
                String[] labelArray = labelList.toArray(new String[0]);

                // CDK grid depiction: 3 columns, auto rows
                dg.withSize(1200, 900) // total image size
                        .depict(molList, 3, -1) // 3 columns, auto rows
                        .writeTo(gridFile.getAbsolutePath());

                System.out.printf("  Grid image saved : %s%n", gridFile.getPath());

            } catch (Exception e) {
                System.out.println("  WARNING: Grid image generation failed: " + e.getMessage());
                System.out.println("  Individual PNG files are still available.");
            }
        }
        // ── Display instructions for NetBeans ─────────────────────────────────

        System.out.printf("  Processed: %d/%d compounds successfully%n",
                successCount, COMPOUNDS.size());

    }
}
