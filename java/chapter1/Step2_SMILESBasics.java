package chemoinformatics;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.exception.InvalidSmilesException;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Step 2 — SMILES Basics
 *
 * SMILES (Simplified Molecular Input Line Entry System) encodes a molecule as a
 * plain text string.
 *
 * SMILES notation rules: - Atoms: C, N, O, H, S, P (element symbols) - Single
 * bond: default (no symbol) - Double bond: = e.g. C=O (formaldehyde) - Triple
 * bond: # e.g. C#N (cyanide) - Branches: () e.g. CC(=O)O (acetic acid) - Rings:
 * digits e.g. c1ccccc1 (benzene) - Charges: [NH4+], [O-] - Aromatic: lowercase
 * c, n, o
 *
 * Reference: Weininger, D. (1988). SMILES, a chemical language and information
 * system. J. Chem. Inf. Comput. Sci., 28(1), 31-36.
 * https://doi.org/10.1021/ci00057a005
 *
 */
public class Step2_SMILESBasics {

    private static final Map<String, String> MOLECULES = new LinkedHashMap<>();

    static {
        MOLECULES.put("Methane", "C");
        MOLECULES.put("Ethanol", "CCO");
        MOLECULES.put("Acetic acid", "CC(=O)O");
        MOLECULES.put("Benzene", "c1ccccc1");
        MOLECULES.put("Nitromethane", "C[N+](=O)[O-]");
        MOLECULES.put("TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        MOLECULES.put("RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        MOLECULES.put("HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        MOLECULES.put("PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        MOLECULES.put("Ammonium nitrate", "[NH4+].[O-][N+]([O-])=O");
        MOLECULES.put("Picric acid", "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        MOLECULES.put("FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
    }

    public static void main(String[] args) {

        // ── Parse and display each molecule ──────────────────────────────────
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // Suppress CDK logging noise to keep console output clean
        System.setProperty("cdk.logging.level", "ERROR");

        System.out.printf("  %-22s  %-10s  %-8s  %s%n",
                "Compound", "Formula", "Atoms", "Parse Status");
        System.out.println("  " + "-".repeat(72));

        int okCount = 0;
        int errorCount = 0;

        for (Map.Entry<String, String> entry : MOLECULES.entrySet()) {
            String name = entry.getKey();
            String smiles = entry.getValue();
            try {
                IAtomContainer mol = parser.parseSmiles(smiles);
                // Configure atom types so implicit H counts are correct
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
                // Get molecular formula
                IMolecularFormula formula
                        = MolecularFormulaManipulator.getMolecularFormula(mol);
                String formulaStr
                        = MolecularFormulaManipulator.getString(formula);
                int atomCount = mol.getAtomCount();
                System.out.printf("  %-22s  %-10s  %-8d  OK%n",
                        name, formulaStr, atomCount);
                okCount++;
            } catch (InvalidSmilesException e) {
                System.out.printf("  %-22s  %-10s  %-8s  ERROR: %s%n",
                        name, "?", "?", e.getMessage());
                errorCount++;
            } catch (Exception e) {
                System.out.printf("  %-22s  %-10s  %-8s  ERROR: %s%n",
                        name, "?", "?", e.getClass().getSimpleName());
                errorCount++;
            }
        }
        // ── Summary ───────────────────────────────────────────────────────────
        System.out.println("  " + "-".repeat(72));
        System.out.printf("  Parsed successfully: %d / %d molecules%n",
                okCount, MOLECULES.size());

    }
}
