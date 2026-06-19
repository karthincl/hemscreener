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
 * SMILES (Simplified Molecular Input Line Entry System) encodes a molecule
 * as a plain text string. This is the universal input format for CDK,
 * just as it is for RDKit in Python.
 *
 * SMILES notation rules:
 *   - Atoms: C, N, O, H, S, P  (element symbols)
 *   - Single bond: default (no symbol)
 *   - Double bond: =   e.g.  C=O  (formaldehyde)
 *   - Triple bond: #   e.g.  C#N  (cyanide)
 *   - Branches: ()    e.g.  CC(=O)O  (acetic acid)
 *   - Rings: digits   e.g.  c1ccccc1  (benzene)
 *   - Charges: [NH4+], [O-]
 *   - Aromatic: lowercase  c, n, o
 *
 * Reference:
 *   Weininger, D. (1988). SMILES, a chemical language and information system.
 *   J. Chem. Inf. Comput. Sci., 28(1), 31-36.
 *   https://doi.org/10.1021/ci00057a005
 *
 * Energetic compound data from:
 *   Klapötke, T.M. (2019). Chemistry of High-Energy Materials, 4th ed.,
 *   Tables 1.1-1.3. De Gruyter. ISBN 978-3-11-062438-6.
 */
public class Step2_SMILESBasics {

    // ── Dataset: compound name -> SMILES ─────────────────────────────────────
    // All compounds are discussed in Klapötke (2019).
    // SMILES strings use the standard Daylight notation.
    private static final Map<String, String> MOLECULES = new LinkedHashMap<>();

    static {
        MOLECULES.put("Methane",           "C");
        MOLECULES.put("Ethanol",           "CCO");
        MOLECULES.put("Acetic acid",       "CC(=O)O");
        MOLECULES.put("Benzene",           "c1ccccc1");
        MOLECULES.put("Nitromethane",      "C[N+](=O)[O-]");
        MOLECULES.put("TNT",               "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        MOLECULES.put("RDX",               "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        MOLECULES.put("HMX",               "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        MOLECULES.put("PETN",              "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        MOLECULES.put("Ammonium nitrate",  "[NH4+].[O-][N+]([O-])=O");
        MOLECULES.put("Picric acid",       "Oc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        MOLECULES.put("FOX-7",             "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
    }

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 2 — SMILES Basics");
        System.out.println("  What SMILES looks like for energetic molecules");
        System.out.println("============================================================\n");

        // Print SMILES notation rules first
        System.out.println("  SMILES NOTATION RULES:");
        System.out.println("  ----------------------");
        System.out.println("  Atoms      : C  N  O  H  S  P  (element symbols)");
        System.out.println("  Double bond: =    e.g.  C=O  (carbonyl)");
        System.out.println("  Triple bond: #    e.g.  C#N  (nitrile)");
        System.out.println("  Branch     : ()   e.g.  CC(=O)O  (acetic acid)");
        System.out.println("  Ring       : num  e.g.  c1ccccc1  (benzene)");
        System.out.println("  Charge     : [NH4+]  [O-]  [N+]");
        System.out.println("  Aromatic   : lowercase  c  n  o");
        System.out.println();

        // ── Parse and display each molecule ──────────────────────────────────
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser       parser  = new SmilesParser(builder);

        // Suppress CDK logging noise to keep console output clean
        System.setProperty("cdk.logging.level", "ERROR");

        System.out.printf("  %-22s  %-10s  %-8s  %s%n",
                          "Compound", "Formula", "Atoms", "Parse Status");
      //  System.out.println("  " + "-".repeat(72));

        int okCount    = 0;
        int errorCount = 0;

        for (Map.Entry<String, String> entry : MOLECULES.entrySet()) {
            String name  = entry.getKey();
            String smiles = entry.getValue();

            try {
                IAtomContainer mol = parser.parseSmiles(smiles);

                // Configure atom types so implicit H counts are correct
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

                // Get molecular formula
                IMolecularFormula formula =
                    MolecularFormulaManipulator.getMolecularFormula(mol);
                String formulaStr =
                    MolecularFormulaManipulator.getString(formula);

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
      //  System.out.println("  " + "-".repeat(72));
        System.out.printf("  Parsed successfully: %d / %d molecules%n",
                          okCount, MOLECULES.size());

        // ── Explain nitro group encoding ──────────────────────────────────────
        System.out.println();
        System.out.println("  NITRO GROUP ENCODING NOTE:");
        System.out.println("  The -NO2 group is written as [N+](=O)[O-] in SMILES.");
        System.out.println("  This zwitterionic form is chemically correct and is the");
        System.out.println("  standard representation used by CDK, RDKit, and PubChem.");
        System.out.println("  Klapötke (2019), Section 1.1 — structural conventions.");
        System.out.println();
        System.out.println("  Proceed to Step 3: Name-to-SMILES lookup via PubChem.");
        System.out.println("============================================================");
    }
}


/*
output: 
run:
============================================================
  STEP 2 — SMILES Basics
  What SMILES looks like for energetic molecules
============================================================

  SMILES NOTATION RULES:
  ----------------------
  Atoms      : C  N  O  H  S  P  (element symbols)
  Double bond: =    e.g.  C=O  (carbonyl)
  Triple bond: #    e.g.  C#N  (nitrile)
  Branch     : ()   e.g.  CC(=O)O  (acetic acid)
  Ring       : num  e.g.  c1ccccc1  (benzene)
  Charge     : [NH4+]  [O-]  [N+]
  Aromatic   : lowercase  c  n  o

  Compound                Formula     Atoms     Parse Status
  Methane                 CH4         5         OK
  Ethanol                 C2H6O       9         OK
  Acetic acid             C2H4O2      8         OK
  Benzene                 C6H6        12        OK
  Nitromethane            CH3NO2      7         OK
  TNT                     C7H5N3O6    21        OK
  RDX                     C3H6N6O6    21        OK
  HMX                     C5H9N8O8    30        OK
  PETN                    C5H8N4O12   29        OK
  Ammonium nitrate        H4N2O3      9         OK
  Picric acid             C6H3N3O7    19        OK
  FOX-7                   C2H4N4O4    14        OK
  Parsed successfully: 12 / 12 molecules

  NITRO GROUP ENCODING NOTE:
  The -NO2 group is written as [N+](=O)[O-] in SMILES.
  This zwitterionic form is chemically correct and is the
  standard representation used by CDK, RDKit, and PubChem.
  Klapötke (2019), Section 1.1 — structural conventions.

  Proceed to Step 3: Name-to-SMILES lookup via PubChem.
============================================================
BUILD SUCCESSFUL (total time: 4 seconds)


*/

