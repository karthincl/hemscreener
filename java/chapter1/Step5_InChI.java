package chemoinformatics;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.interfaces.IMolecularFormula;

import java.util.LinkedHashMap;
import java.util.Map;

public class Step5_InChI {

    public static class MoleculeIdentifiers {

        public final String name;
        public final String smiles;
        public final String canonicalSmiles;
        public final String formula;
        public final double molecularWeight;
        public final String inchi;
        public final String inchiKey;
        public final String inchiStatus;   // "OK", "WARNING", "ERROR"

        public MoleculeIdentifiers(String name, String smiles, String canonical,
                String formula, double mw,
                String inchi, String inchiKey, String status) {
            this.name = name;
            this.smiles = smiles;
            this.canonicalSmiles = canonical;
            this.formula = formula;
            this.molecularWeight = mw;
            this.inchi = inchi;
            this.inchiKey = inchiKey;
            this.inchiStatus = status;
        }
    }

    // ── Core method: compute all identifiers from SMILES ─────────────────────
    public static MoleculeIdentifiers computeIdentifiers(
            String name, String smiles,
            SmilesParser parser,
            SmilesGenerator smiGen,
            InChIGeneratorFactory inchiFactory) {

        try {
            // 1. Parse SMILES
            IAtomContainer mol = parser.parseSmiles(smiles);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
            // 2. Canonical SMILES
            String canonical = smiGen.create(mol);
            // 3. Molecular formula and weight
            IMolecularFormula mfObj
                    = MolecularFormulaManipulator.getMolecularFormula(mol);
            String formulaStr
                    = MolecularFormulaManipulator.getString(mfObj);
            double mw
                    = MolecularFormulaManipulator.getTotalExactMass(mfObj);
            // 4. InChI generation via CDK's JNA-InChI bridge
            InChIGenerator gen = inchiFactory.getInChIGenerator(mol);
            String inchiStr = gen.getInchi();
            String inchiKeyStr = gen.getInchiKey();
            String returnCode = gen.getReturnStatus().toString();
            return new MoleculeIdentifiers(name, smiles, canonical,
                    formulaStr, mw,
                    inchiStr, inchiKeyStr, returnCode);
        } catch (Exception e) {
            return new MoleculeIdentifiers(name, smiles, "ERROR", "?", 0.0,
                    "ERROR: " + e.getMessage(), "?", "ERROR");
        }
    }

    public static void main(String[] args) {
        // ── Setup CDK objects ─────────────────────────────────────────────────
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);
        SmilesGenerator smiGen = SmilesGenerator.unique();
        InChIGeneratorFactory inchiFactory;
        try {
            inchiFactory = InChIGeneratorFactory.getInstance();
        } catch (Exception e) {
            System.err.println("ERROR: Could not initialize InChI factory.");
            System.err.println("  Make sure cdk-2.9.jar is on the classpath.");
            System.err.println("  Details: " + e.getMessage());
            return;
        }
        // ── Dataset: key HEMs from Klapötke (2019) ───────────────────────────
        Map<String, String> compounds = new LinkedHashMap<>();
        compounds.put("Nitromethane", "C[N+](=O)[O-]");
        compounds.put("TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        compounds.put("RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        compounds.put("HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        compounds.put("PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        compounds.put("Picric acid", "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        compounds.put("FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        compounds.put("NTO", "O=c1[nH]nno1");
        compounds.put("TATB", "Nc1c([N+](=O)[O-])nc([N+](=O)[O-])nc1[N+](=O)[O-]");
        // ── Compute and display identifiers ───────────────────────────────────
        System.out.println("  Full identifier set for each compound:");
        System.out.println("  " + "=".repeat(60));

        for (Map.Entry<String, String> entry : compounds.entrySet()) {
            MoleculeIdentifiers ids
                    = computeIdentifiers(entry.getKey(), entry.getValue(),
                            parser, smiGen, inchiFactory);
            System.out.println();
            System.out.println("  " + ids.name);
            System.out.println("  " + "-".repeat(50));
            System.out.printf("    Formula         : %s%n", ids.formula);
            System.out.printf("    Exact MW        : %.5f g/mol%n", ids.molecularWeight);
            System.out.printf("    Input SMILES    : %s%n", ids.smiles);
            System.out.printf("    Canonical SMILES: %s%n", ids.canonicalSmiles);

            // Display InChI with line break after the formula layer for readability
            if (ids.inchi != null && ids.inchi.startsWith("InChI=")) {
                System.out.printf("    InChI           : %s%n", ids.inchi);
            } else {
                System.out.printf("    InChI           : %s%n", ids.inchi);
            }
            System.out.printf("    InChIKey        : %s%n", ids.inchiKey);
            System.out.printf("    InChI status    : %s%n", ids.inchiStatus);
        }
        
    }
}
