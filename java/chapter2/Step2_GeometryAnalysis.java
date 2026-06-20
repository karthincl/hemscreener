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
//  Chapter 2 — Step 2: Bond Lengths and Angles Analysis
//  Chemoinformatics & Computational Modeling of Organic Molecules
//  Reference:
//  M.Karthikeyan, Renu Vyas (2014). Practical Chemoinformatics, Springer.
// ============================================================
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.modeling.builder3d.TemplateHandler3D;
import org.openscience.cdk.geometry.GeometryUtil;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Chapter 2, Step 2 — Bond Lengths and Angles Analysis.
 *
 * BACKGROUND: ----------- After generating 3D coordinates (Step 1), the first
 * validation check is to inspect bond lengths and angles against known
 * reference values. Unrealistic values indicate a failed geometry generation.
 *
 * KEY BOND LENGTH REFERENCES FOR ENERGETIC MOLECULES: (from Klapötke 2019,
 * Appendix A; Pauling 1960)
 *
 * Bond Typical length (Å) Notes C-C single 1.54 alkane C=C double 1.34 alkene
 * C-N single 1.47 amine N-N single 1.45 hydrazine-type N=N double 1.25 azo N-O
 * single 1.36 hydroxylamine N=O double 1.22 nitro (N=O) N-NO2 1.39 nitramine
 * (weakest, key sensitivity bond) O-NO2 1.40 nitrate ester C-NO2 1.47
 * nitroalkane
 *
 * The N-NO2 bond (nitramine bond) is the weakest bond in RDX and HMX. Its bond
 * dissociation energy (~147 kJ/mol, CBS-QB3) is the primary sensitivity
 * indicator. Klapötke (2019), Chapter 3.
 *
 * BOND ANGLE REFERENCES: sp3 carbon 109.5° sp2 carbon 120° N in nitro ~117°
 * (O-N-O), ~121° (C-N=O or N-N=O)
 *
 * REFERENCES: ----------- 1. Klapötke, T.M. (2019). Chemistry of High-Energy
 * Materials, 4th ed. Appendix A (bond lengths), Chapter 3 (sensitivity). De
 * Gruyter. 2. Pauling, L. (1960). The Nature of the Chemical Bond, 3rd ed.
 * Cornell University Press.
 */
public class Step2_GeometryAnalysis {

    // ── Reference bond lengths (Angstrom) from Klapötke (2019) Appendix ─────
    // Used to flag unrealistic values in generated geometries
    private static final double CC_SINGLE_REF = 1.54;
    private static final double CN_SINGLE_REF = 1.47;
    private static final double NN_SINGLE_REF = 1.45;
    private static final double NO_DOUBLE_REF = 1.22;   // N=O in nitro group
    private static final double NNO2_REF = 1.39;   // nitramine N-NO2 bond
    private static final double TOLERANCE = 0.15;   // +/- 0.15 Å acceptable

    // ── Result data classes ───────────────────────────────────────────────────
    public static class BondResult {

        public final String atom1Symbol;
        public final int atom1Idx;
        public final String atom2Symbol;
        public final int atom2Idx;
        public final double bondOrder;
        public final double length;       // Angstrom
        public final String bondType;     // descriptive label
        public final boolean flagged;     // true if length seems unrealistic

        public BondResult(String s1, int i1, String s2, int i2,
                double order, double len, String type, boolean flag) {
            this.atom1Symbol = s1;
            this.atom1Idx = i1;
            this.atom2Symbol = s2;
            this.atom2Idx = i2;
            this.bondOrder = order;
            this.length = len;
            this.bondType = type;
            this.flagged = flag;
        }
    }

    public static class AngleResult {

        public final String atom1Symbol;
        public final String centralSymbol;
        public final String atom3Symbol;
        public final int atom1Idx;
        public final int centralIdx;
        public final int atom3Idx;
        public final double angleDeg;

        public AngleResult(String s1, int i1, String sc, int ic,
                String s3, int i3, double angle) {
            this.atom1Symbol = s1;
            this.atom1Idx = i1;
            this.centralSymbol = sc;
            this.centralIdx = ic;
            this.atom3Symbol = s3;
            this.atom3Idx = i3;
            this.angleDeg = angle;
        }
    }

    // ── Bond length calculation ───────────────────────────────────────────────
    /**
     * Calculate the distance between two atoms in 3D space. Equivalent to
     * rdMolTransforms.GetBondLength() in Python RDKit.
     *
     * @param mol molecule with 3D coordinates
     * @param idx1 index of first atom
     * @param idx2 index of second atom
     * @return distance in Angstrom
     */
    public static double bondLength(IAtomContainer mol, int idx1, int idx2) {
        Point3d p1 = mol.getAtom(idx1).getPoint3d();
        Point3d p2 = mol.getAtom(idx2).getPoint3d();
        if (p1 == null || p2 == null) {
            return -1.0;
        }
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    // ── Bond angle calculation ────────────────────────────────────────────────
    /**
     * Calculate the bond angle at the central atom (in degrees). Equivalent to
     * rdMolTransforms.GetAngleDeg() in Python RDKit.
     *
     * Computes the angle between vectors (central->atom1) and (central->atom3).
     *
     * @param mol molecule with 3D coordinates
     * @param idx1 index of first atom
     * @param idxC index of central atom
     * @param idx3 index of third atom
     * @return angle in degrees
     */
    public static double bondAngleDeg(IAtomContainer mol,
            int idx1, int idxC, int idx3) {
        Point3d p1 = mol.getAtom(idx1).getPoint3d();
        Point3d pC = mol.getAtom(idxC).getPoint3d();
        Point3d p3 = mol.getAtom(idx3).getPoint3d();
        if (p1 == null || pC == null || p3 == null) {
            return -1.0;
        }

        // Vectors from central atom to the two neighbours
        double[] v1 = {p1.x - pC.x, p1.y - pC.y, p1.z - pC.z};
        double[] v3 = {p3.x - pC.x, p3.y - pC.y, p3.z - pC.z};

        double dot = v1[0] * v3[0] + v1[1] * v3[1] + v1[2] * v3[2];
        double mag1 = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        double mag3 = Math.sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);

        if (mag1 < 1e-10 || mag3 < 1e-10) {
            return -1.0;
        }

        double cosAngle = dot / (mag1 * mag3);
        // Clamp to [-1, 1] to avoid NaN from floating-point errors
        cosAngle = Math.max(-1.0, Math.min(1.0, cosAngle));
        return Math.toDegrees(Math.acos(cosAngle));
    }

    // ── Analyse all bonds in a molecule ──────────────────────────────────────
    /**
     * Extract and classify all non-H bonds with their lengths.
     *
     * @param mol molecule with 3D coordinates
     * @return list of BondResult objects, one per non-H bond
     */
    public static List<BondResult> analyseBonds(IAtomContainer mol) {
        List<BondResult> results = new ArrayList<BondResult>();

        for (IBond bond : mol.bonds()) {
            IAtom a1 = bond.getAtom(0);
            IAtom a2 = bond.getAtom(1);

            // Skip bonds involving hydrogen for cleaner output
            if (a1.getSymbol().equals("H") || a2.getSymbol().equals("H")) {
                continue;
            }

            int i1 = mol.indexOf(a1);
            int i2 = mol.indexOf(a2);
            double len = bondLength(mol, i1, i2);
            double order = bond.getOrder() != null
                    ? bond.getOrder().numeric() : 1.0;

            String s1 = a1.getSymbol();
            String s2 = a2.getSymbol();

            // Classify bond type with sensitivity annotation
            String btype = classifyBond(s1, s2, order, bond, mol, i1, i2);

            // Flag unrealistic lengths
            boolean flagged = (len > 0) && (len < 0.8 || len > 2.2);

            results.add(new BondResult(s1, i1, s2, i2, order, len, btype, flagged));
        }
        return results;
    }

    // ── Bond type classifier ──────────────────────────────────────────────────
    private static String classifyBond(String s1, String s2, double order,
            IBond bond, IAtomContainer mol,
            int i1, int i2) {
        // Nitramine N-NO2 bond: N bonded to another N that has two O neighbours
        if (s1.equals("N") && s2.equals("N")) {
            if (isNitroNitrogen(mol, i2) || isNitroNitrogen(mol, i1)) {
                return "N-NO2 nitramine [SENSITIVITY KEY]";
            }
            if (order >= 2.0) {
                return "N=N azo";
            }
            return "N-N single";
        }
        // Nitro group N-O bonds
        if ((s1.equals("N") && s2.equals("O"))
                || (s1.equals("O") && s2.equals("N"))) {
            if (order >= 2.0) {
                return "N=O nitro";
            }
            return "N-O single";
        }
        // Carbon-nitrogen
        if ((s1.equals("C") && s2.equals("N"))
                || (s1.equals("N") && s2.equals("C"))) {
            if (mol.getAtom(s1.equals("N") ? i1 : i2).isAromatic()) {
                return "C-N aromatic";
            }
            if (order >= 2.0) {
                return "C=N imine";
            }
            return "C-N single";
        }
        // Carbon-carbon
        if (s1.equals("C") && s2.equals("C")) {
            if (order >= 2.0) {
                return "C=C double";
            }
            if (mol.getAtom(i1).isAromatic()) {
                return "C-C aromatic";
            }
            return "C-C single";
        }
        // Carbon-oxygen
        if ((s1.equals("C") && s2.equals("O"))
                || (s1.equals("O") && s2.equals("C"))) {
            if (order >= 2.0) {
                return "C=O carbonyl";
            }
            return "C-O single";
        }
        // Oxygen-nitrogen nitrate ester O-NO2
        if ((s1.equals("O") && s2.equals("N"))
                || (s1.equals("N") && s2.equals("O"))) {
            return "O-N nitrate ester";
        }
        return s1 + "-" + s2;
    }

    /**
     * Returns true if atom at idx is a nitrogen with >= 2 oxygen neighbours.
     */
    private static boolean isNitroNitrogen(IAtomContainer mol, int idx) {
        IAtom n = mol.getAtom(idx);
        if (!n.getSymbol().equals("N")) {
            return false;
        }
        int oCount = 0;
        for (IBond b : mol.getConnectedBondsList(n)) {
            IAtom other = b.getOther(n);
            if (other.getSymbol().equals("O")) {
                oCount++;
            }
        }
        return oCount >= 2;
    }

    // ── Analyse all nitrogen-centred angles ───────────────────────────────────
    /**
     * Calculate bond angles around nitrogen atoms only. Nitrogen geometry is
     * diagnostic for the type of N-functional group: sp3 N (amine) -> ~109.5°
     * sp2 N (nitro, imine) -> ~120° planar N-NO2 -> ~117° (O-N-O)
     */
    public static List<AngleResult> analyseNitrogenAngles(IAtomContainer mol) {
        List<AngleResult> results = new ArrayList<AngleResult>();

        for (int c = 0; c < mol.getAtomCount(); c++) {
            IAtom central = mol.getAtom(c);
            if (!central.getSymbol().equals("N")) {
                continue;
            }

            // Get non-H neighbours
            List<Integer> nbrs = new ArrayList<Integer>();
            for (IBond b : mol.getConnectedBondsList(central)) {
                IAtom other = b.getOther(central);
                if (!other.getSymbol().equals("H")) {
                    nbrs.add(mol.indexOf(other));
                }
            }

            // All pairs of neighbours
            for (int i = 0; i < nbrs.size(); i++) {
                for (int j = i + 1; j < nbrs.size(); j++) {
                    double angle = bondAngleDeg(mol, nbrs.get(i), c, nbrs.get(j));
                    if (angle > 0) {
                        results.add(new AngleResult(
                                mol.getAtom(nbrs.get(i)).getSymbol(), nbrs.get(i),
                                "N", c,
                                mol.getAtom(nbrs.get(j)).getSymbol(), nbrs.get(j),
                                angle));
                    }
                }
            }
        }
        return results;
    }

    // ── Main ─────────────────────────────────────────────────────────────────
    public static void main(String[] args) {

        System.out.println(rep('=', 62));
        System.out.println("  CHAPTER 2 - STEP 2: Bond Lengths and Angles Analysis");
        System.out.println("  Based on: Klapötke (2019), Chemistry of High-Energy");
        System.out.println("  Materials, 4th ed., Appendix A & Chapter 3");
        System.out.println(rep('=', 62));
        System.out.println();

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // ── Reference table of bond lengths ──────────────────────────────────
        System.out.println("  REFERENCE BOND LENGTHS (Klapötke 2019, Appendix A):");
        System.out.println("  " + rep('-', 50));
        System.out.printf("  %-28s  %s%n", "Bond Type", "Typical Length (Å)");
        System.out.println("  " + rep('-', 50));
        String[][] refs = {
            {"C-C single (alkane)", "1.54"},
            {"C=C double (alkene)", "1.34"},
            {"C-N single", "1.47"},
            {"N-N single", "1.45"},
            {"N=N double (azo)", "1.25"},
            {"N=O (nitro group)", "1.22"},
            {"N-NO2 (nitramine)", "1.39  <-- weakest bond in RDX/HMX"},
            {"O-NO2 (nitrate ester)", "1.40"},
            {"C-NO2 (nitroalkane)", "1.47"},};
        for (String[] r : refs) {
            System.out.printf("  %-28s  %s%n", r[0], r[1]);
        }
        System.out.println();

        // ── Analyse TNT ───────────────────────────────────────────────────────
        analyseCompound("TNT",
                "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",
                parser, true, true);

        // ── Analyse RDX ───────────────────────────────────────────────────────
        analyseCompound("RDX",
                "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]",
                parser, true, true);

        // ── Analyse PETN ──────────────────────────────────────────────────────
        analyseCompound("PETN",
                "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]",
                parser, true, false);

        System.out.println("  KEY INTERPRETATION:");
        System.out.println("  - N-NO2 bonds in RDX/HMX are the SHORTEST N-N bonds.");
        System.out.println("    Shorter bond order -> weaker homolytic BDE (~147 kJ/mol).");
        System.out.println("  - C-NO2 bonds in TNT are LONGER and STRONGER (~258 kJ/mol).");
        System.out.println("    This is why TNT is LESS sensitive than RDX.");
        System.out.println("  Reference: Klapötke (2019), Table 3.1 — BDE values.");
        System.out.println(rep('=', 62));
    }

    // ── Per-compound analysis routine ────────────────────────────────────────
    private static void analyseCompound(String name, String smiles,
            SmilesParser parser,
            boolean showBonds,
            boolean showAngles) {
        System.out.println("  " + rep('-', 62));
        System.out.println("  COMPOUND: " + name);
        System.out.println("  SMILES  : " + smiles);
        System.out.println("  " + rep('-', 62));

        // Generate 3D structure
        Step1_3DGeneration.Structure3D s
                = Step1_3DGeneration.generateStructure3D(name, smiles, parser);

        if (!s.success || s.mol == null) {
            System.out.println("  3D generation failed: " + s.message);
            System.out.println();
            return;
        }
        System.out.println("  3D status: " + s.message);

        // ── Bond lengths ──────────────────────────────────────────────────────
        if (showBonds) {
            List<BondResult> bonds = analyseBonds(s.mol);
            System.out.println();
            System.out.printf("  %-6s  %-6s  %-8s  %-12s  %s%n",
                    "Atom1", "Atom2", "Order", "Length (Å)", "Bond Type");
            System.out.println("  " + rep('-', 65));

            for (BondResult b : bonds) {
                String flagStr = b.flagged ? " [CHECK]" : "";
                System.out.printf("  %s%-3d  %s%-3d  %-8.1f  %-12.4f  %s%s%n",
                        b.atom1Symbol, b.atom1Idx,
                        b.atom2Symbol, b.atom2Idx,
                        b.bondOrder, b.length,
                        b.bondType, flagStr);
            }
        }

        // ── Nitrogen angles ───────────────────────────────────────────────────
        if (showAngles) {
            List<AngleResult> angles = analyseNitrogenAngles(s.mol);
            System.out.println();
            System.out.printf("  Bond Angles around Nitrogen atoms:%n");
            System.out.printf("  %-18s  %s%n", "Angle", "Degrees");
            System.out.println("  " + rep('-', 35));

            for (AngleResult a : angles) {
                String label = String.format("%s%d-N%d-%s%d",
                        a.atom1Symbol, a.atom1Idx,
                        a.centralIdx,
                        a.atom3Symbol, a.atom3Idx);
                System.out.printf("  %-18s  %.2f%n", label, a.angleDeg);
            }
        }
        System.out.println();
    }

    static String rep(char c, int n) {
        return new String(new char[n]).replace('\0', c);
    }
}
