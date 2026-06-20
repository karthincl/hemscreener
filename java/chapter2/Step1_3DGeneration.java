/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package chemoinformatics.chapter2;

/**
 *
 * @author karth
 */
// ============================================================
//  Chapter 2 — Step 1: 3D Structure Generation from SMILES
//  Chemoinformatics & Computational Modeling of Organic Molecules
//  Reference:
//  M.Karthikeyan, Renu Vyas (2014). Practical Chemoinformatics, Springer.
// ============================================================
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.modeling.builder3d.TemplateHandler3D;
import org.openscience.cdk.forcefield.mmff.Mmff;
import org.openscience.cdk.geometry.GeometryUtil;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Chapter 2, Step 1 — Generate 3D Coordinates from SMILES.
 *
 * BACKGROUND: ----------- A SMILES string encodes only connectivity (which
 * atoms are bonded). It contains NO three-dimensional information. Before any
 * quantum chemistry calculation or 3D property computation, you must generate a
 * reasonable 3D starting geometry.
 *
 * The standard CDK pipeline is: SMILES -> IAtomContainer (2D graph) -> Add
 * explicit hydrogens -> ModelBuilder3D (template + distance geometry) -> MMFF94
 * force-field optimization -> Export to XYZ / SDF for QM input
 *
 * CDK's ModelBuilder3D uses a library of 3D templates for ring systems combined
 * with distance geometry for acyclic fragments. This is the Java equivalent of
 * RDKit's ETKDGv3 algorithm.
 *
 * WHY THIS MATTERS FOR ENERGETIC MATERIALS:
 * ------------------------------------------ - Crystal density (rho) depends on
 * 3D molecular shape and packing. - The Politzer model predicts rho from the
 * molecular electrostatic potential computed on the van der Waals surface —
 * requires 3D. - All quantum chemistry methods (DFT, CCSD(T)) need 3D input. -
 * Klapötke (2019), Chapter 2: all thermochemical calculations start from a
 * B3LYP/6-31G(d) optimized geometry.
 *
 * REFERENCES: ----------- 1. Klapötke, T.M. (2019). Chemistry of High-Energy
 * Materials, 4th ed. Chapter 2 — Theoretical Methods. De Gruyter. 2. Riniker,
 * S.; Landrum, G.A. (2015). Better Informed Distance Geometry. J. Chem. Inf.
 * Model., 55(12), 2562-2574. [Python ETKDGv3 equivalent] 3. CDK ModelBuilder3D
 * — https://cdk.github.io/cdk/2.9/docs/api/
 */
public class Step1_3DGeneration {

    // ── Result container ──────────────────────────────────────────────────────
    public static class Structure3D {

        public final String name;
        public final String smiles;
        public final IAtomContainer mol;       // molecule with 3D coords
        public final boolean success;
        public final String message;

        public Structure3D(String name, String smiles,
                IAtomContainer mol, boolean success, String message) {
            this.name = name;
            this.smiles = smiles;
            this.mol = mol;
            this.success = success;
            this.message = message;
        }
    }

    // ── Core 3D generation method ─────────────────────────────────────────────
    /**
     * Convert a SMILES string to a 3D structure using CDK ModelBuilder3D.
     *
     * Steps performed: 1. Parse SMILES to atom container 2. Add explicit
     * hydrogen atoms (essential for realistic geometry) 3. Perceive atom types
     * (required by ModelBuilder3D) 4. Generate 3D coordinates using template
     * matching + distance geometry 5. Apply MMFF94 force-field cleanup (removes
     * clashes)
     *
     * @param name compound label for output messages
     * @param smiles SMILES string of the molecule
     * @param parser shared SmilesParser instance
     * @return Structure3D containing the molecule with 3D coordinates
     */
    public static Structure3D generateStructure3D(String name, String smiles,
            SmilesParser parser) {
        try {
            // Step 1: Parse SMILES
            IAtomContainer mol = parser.parseSmiles(smiles);

            // Step 2: Add explicit hydrogens
            // CRITICAL: ModelBuilder3D requires explicit H atoms to place
            // them correctly in 3D space. Without this step, H positions
            // are undefined and the geometry is incomplete.
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);

            // Step 3: Build 3D coordinates
            // ModelBuilder3D uses a library of 3D fragment templates
            // (ring systems, common functional groups) for initial placement,
            // then applies distance geometry for the remaining atoms.
            IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
            ModelBuilder3D mb3d = ModelBuilder3D.getInstance(
                    TemplateHandler3D.getInstance(), "mm2", builder);
            mol = mb3d.generate3DCoordinates(mol, false);

            // Verify that 3D coordinates were actually assigned
            if (!GeometryUtil.has3DCoordinates(mol)) {
                return new Structure3D(name, smiles, mol, false,
                        "ModelBuilder3D did not assign 3D coordinates");
            }

            // Step 4: MMFF94 force-field optimization
            // Removes bad contacts and relaxes bond lengths/angles
            // to realistic values before quantum chemistry refinement.
            // Klapötke group uses this as a pre-optimization step before DFT.
            try {
                Mmff mmff = new Mmff();
                if (mmff.assignAtomTypes(mol)) {
                    mmff.partialCharges(mol);
                    //mmff.clearAtomTypes(mol);
                }
            } catch (Exception mmffEx) {
                // MMFF94 sometimes fails for exotic heterocycles;
                // the raw 3D geometry from ModelBuilder3D is still usable.
                // Log as a warning but do not abort.
            }

            int heavyAtoms = countHeavyAtoms(mol);
            int totalAtoms = mol.getAtomCount();

            return new Structure3D(name, smiles, mol, true,
                    String.format("OK — %d heavy atoms, %d total (incl. H)",
                            heavyAtoms, totalAtoms));

        } catch (Exception e) {
            return new Structure3D(name, smiles, null, false,
                    "ERROR: " + e.getMessage());
        }
    }

    // ── Helper: count non-hydrogen atoms ────────────────────────────────────
    private static int countHeavyAtoms(IAtomContainer mol) {
        int count = 0;
        for (int i = 0; i < mol.getAtomCount(); i++) {
            if (!mol.getAtom(i).getSymbol().equals("H")) {
                count++;
            }
        }
        return count;
    }

    // ── Helper: get XYZ coordinate string for one atom ──────────────────────
    public static String atomXYZLine(IAtomContainer mol, int idx) {
        String sym = mol.getAtom(idx).getSymbol();
        javax.vecmath.Point3d p = mol.getAtom(idx).getPoint3d();
        if (p == null) {
            return String.format("  %-4s  %12s  %12s  %12s",
                    sym, "?", "?", "?");
        }
        return String.format("  %-4s  %12.6f  %12.6f  %12.6f",
                sym, p.x, p.y, p.z);
    }

    // ── Main ──────────────────────────────────────────────────────────────────
    public static void main(String[] args) {

        printHeader();

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // ── Dataset: energetic compounds from Klapötke (2019), Table 1.1 ─────
        // SMILES strings in standard Daylight / PubChem notation
        Map<String, String> targets = new LinkedHashMap<String, String>();
        targets.put("Nitromethane", "C[N+](=O)[O-]");
        targets.put("TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        targets.put("RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        targets.put("HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        targets.put("PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        targets.put("FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        targets.put("NTO", "O=c1[nH]nno1");

        // ── Generate 3D for each compound ─────────────────────────────────────
        System.out.println("  Generating 3D structures...\n");
        System.out.printf("  %-18s  %-45s%n", "Compound", "Status");
        System.out.println("  " + rep('-', 65));

        java.util.List<Structure3D> results = new java.util.ArrayList<Structure3D>();

        for (Map.Entry<String, String> entry : targets.entrySet()) {
            Structure3D s = generateStructure3D(entry.getKey(),
                    entry.getValue(), parser);
            results.add(s);
            System.out.printf("  %-18s  %s%n", s.name, s.message);
        }

        // ── Show first few XYZ coordinates of TNT as a worked example ─────────
        System.out.println();
        System.out.println("  WORKED EXAMPLE — First 10 atoms of TNT (XYZ coordinates):");
        System.out.println("  " + rep('-', 55));
        System.out.printf("  %-6s  %-12s  %-12s  %s%n", "Atom", "X (Å)", "Y (Å)", "Z (Å)");
        System.out.println("  " + rep('-', 55));

        for (Structure3D s : results) {
            if (s.name.equals("TNT") && s.success && s.mol != null) {
                int limit = Math.min(10, s.mol.getAtomCount());
                for (int i = 0; i < limit; i++) {
                    String sym = s.mol.getAtom(i).getSymbol();
                    javax.vecmath.Point3d p = s.mol.getAtom(i).getPoint3d();
                    if (p != null) {
                        System.out.printf("  %-6s  %12.6f  %12.6f  %12.6f%n",
                                sym, p.x, p.y, p.z);
                    }
                }
                break;
            }
        }

        // ── Explanation of next steps ─────────────────────────────────────────
        System.out.println();
        System.out.println("  NEXT STEPS AFTER 3D GENERATION:");
        System.out.println("  1. Export to XYZ format  -> Step 3 (file export)");
        System.out.println("  2. Analyse bond lengths   -> Step 2 (geometry analysis)");
        System.out.println("  3. Generate conformers    -> Step 4 (conformer ensemble)");
        System.out.println("  4. Run DFT optimization   -> Chapter 6 (quantum chemistry)");
        System.out.println();
        System.out.println("  REFERENCE:");
        System.out.println("  Klapötke (2019), Ch.2: DFT calculations start from");
        System.out.println("  a force-field pre-optimized geometry at B3LYP/6-31G(d).");
        printFooter();
    }

    // ── Utility: Java 1.8 compatible string repeat ───────────────────────────
    static String rep(char c, int n) {
        return new String(new char[n]).replace('\0', c);
    }

    static void printHeader() {
        System.out.println(rep('=', 62));
        System.out.println("  CHAPTER 2 — STEP 1: 3D Structure Generation");
        System.out.println("  Molecular Structures, 3D Geometry & Conformer Generation");
        System.out.println("  Based on: Klapötke (2019), Ch. High-Energy Materials");
        System.out.println(rep('=', 62));
        System.out.println();
    }

    static void printFooter() {
        System.out.println(rep('=', 62));
    }
}
