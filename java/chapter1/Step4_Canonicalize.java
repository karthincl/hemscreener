package chemoinformatics;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.exception.InvalidSmilesException;

public class Step4_Canonicalize {

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 4 — Canonical SMILES Normalization");
        System.out.println("============================================================\n");

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // CDK's canonical SMILES generator (unique() applies Morgan algorithm)
        SmilesGenerator smiGen = SmilesGenerator.unique();

        // ── Part A: Show that different SMILES give the same canonical form ───
        System.out.println("  PART A: Different inputs -> same canonical SMILES\n");
        System.out.println("  Compound: Acetic acid");
        System.out.println();

        String[] aceticAcidVariants = {
            "CC(=O)O", // standard textbook notation
            "OC(C)=O", // reverse traversal
            "C(C)(=O)O", // explicit branch
            "CC(O)=O", // O before =O
        };

        System.out.printf("  %-28s  %-30s  %s%n",
                "Input SMILES", "Canonical SMILES", "Valid?");
        System.out.println("  " + "-".repeat(70));

        for (String smi : aceticAcidVariants) {
            try {
                IAtomContainer mol = parser.parseSmiles(smi);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                String canonical = smiGen.create(mol);
                System.out.printf("  %-28s  %-30s  YES%n", smi, canonical);
            } catch (Exception e) {
                System.out.printf("  %-28s  %-30s  NO  (%s)%n",
                        smi, "(parse failed)", e.getMessage());
            }
        }

        // ── Part B: Normalize SMILES for key energetic compounds ──────────────
        System.out.println("\n  PART B: Canonical SMILES for energetic compounds");
        System.out.println("  Source: Klapötke (2019), Table 1.1\n");

        // Each entry: {name, raw SMILES from PubChem or literature}
        // Some have non-canonical forms to demonstrate normalization
        String[][] hemCompounds = {
            {"Nitromethane", "C[N+](=O)[O-]"},
            {"TNT", "[N+](=O)([O-])c1cc([N+](=O)[O-])cc(c1C)[N+](=O)[O-]"},
            {"RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]"},
            {"HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1"},
            {"PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]"},
            {"Picric acid", "c1c(O)c([N+](=O)[O-])cc([N+](=O)[O-])c1[N+](=O)[O-]"},
            {"FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N"},
            {"NTO", "O=c1[nH]nno1"},
            {"TATB", "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]"},
            {"Ammonium nitrate", "[NH4+].[O-][N+]([O-])=O"},};

        System.out.printf("  %-22s  %-14s  %s%n",
                "Name", "Formula", "Canonical SMILES");
        System.out.println("  " + "-".repeat(80));

        for (String[] entry : hemCompounds) {
            String name = entry[0];
            String raw = entry[1];
            try {
                IAtomContainer mol = parser.parseSmiles(raw);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

                IMolecularFormula formula
                        = MolecularFormulaManipulator.getMolecularFormula(mol);
                String formulaStr
                        = MolecularFormulaManipulator.getString(formula);
                String canonical = smiGen.create(mol);
                // Truncate long SMILES for display
                String display = canonical.length() > 42
                        ? canonical.substring(0, 42) + "..."
                        : canonical;
                System.out.printf("  %-22s  %-14s  %s%n", name, formulaStr, display);
            } catch (Exception e) {
                System.out.printf("  %-22s  %-14s  ERROR: %s%n",
                        name, "?", e.getMessage());
            }
        }
        // ── Part C: Demonstrate stereo-specific SMILES ────────────────────────
      
        SmilesGenerator isoGen = SmilesGenerator.isomeric();
        String[] stereoSmiles = {
            "C[C@@H](O)N", // L-alaninol — stereo specified
            "C[C@H](O)N", // D-alaninol — opposite stereo
        };

        System.out.printf("  %-22s  %-22s  %s%n",
                "Input", "Canonical (no stereo)", "Isomeric (stereo)");
        System.out.println("  " + "-".repeat(70));

        for (String smi : stereoSmiles) {
            try {
                IAtomContainer mol = parser.parseSmiles(smi);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                String canonical = smiGen.create(mol);
                String isomeric = isoGen.create(mol);
                System.out.printf("  %-22s  %-22s  %s%n", smi, canonical, isomeric);
            } catch (Exception e) {
                System.out.printf("  %-22s  ERROR: %s%n", smi, e.getMessage());
            }
        }

    }
}
