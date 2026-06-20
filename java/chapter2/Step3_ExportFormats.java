/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
/**
 *
 * @author karth
 */
package chemoinformatics.chapter2;

// ============================================================
//  Chapter 2 — Step 3: Export 3D Structures to File Formats
//  Chemoinformatics & Computational Modeling of Organic Molecules
//  Chemoinformatics & Computational Modeling of Organic Molecules
//  Reference:
//  M.Karthikeyan, Renu Vyas (2014). Practical Chemoinformatics, Springer.
//  Output files created in: structures/ subfolder
//    TNT.sdf, TNT.mol, TNT.xyz
//    RDX.sdf, RDX.mol, RDX.xyz
//    ...etc
// ============================================================
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.SDFWriter;

import javax.vecmath.Point3d;
import java.io.File;
import java.io.FileWriter;
import java.io.StringWriter;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Chapter 2, Step 3 — Export 3D Structures to SDF, XYZ and MOL Formats.
 *
 * BACKGROUND — File Formats in Computational Chemistry:
 * -------------------------------------------------------
 *
 * SDF (Structure-Data File, .sdf): The universal interchange format. Supports
 * multiple molecules per file (separated by $$$$). Contains atom coordinates,
 * bond table, and optional data fields. Used by all cheminformatics tools
 * (RDKit, CDK, OpenBabel, MOE, Schrodinger). CDK API: SDFWriter
 *
 * MOL (MDL Molfile, .mol): Single-structure SDF. The V2000 format is supported
 * by all tools. CDK API: MDLV2000Writer
 *
 * XYZ (Cartesian coordinate file, .xyz): The input format for ORCA and Gaussian
 * quantum chemistry programs. Klapötke's group starts all DFT calculations from
 * XYZ files pre-optimized with MMFF94 (Klapötke 2019, Chapter 2). Format:
 * line1=atom count, line2=title, then one "Symbol X Y Z" per atom. CDK API:
 * write manually (simple format, no library needed)
 *
 * WORKFLOW IN KLAPÖTKE'S GROUP: SMILES -> CDK 3D -> XYZ -> ORCA B3LYP/6-31G(d)
 * -> CBS-4M -> ΔHf
 *
 * REFERENCES: ----------- 1. Klapötke, T.M. (2019). Chemistry of High-Energy
 * Materials, 4th ed. Chapter 2 — Theoretical Methods. De Gruyter. 2. Dalby, A.
 * et al. (1992). Description of several chemical structure file formats used by
 * computer programs developed at Molecular Design Limited. J. Chem. Inf.
 * Comput. Sci., 32(3), 244-255. 3. Neese, F. et al. (2020). The ORCA program
 * system. J. Chem. Phys., 152, 224108. (XYZ input format)
 */
public class Step3_ExportFormats {

    // ── Export result container ───────────────────────────────────────────────
    public static class ExportResult {

        public final String name;
        public final String sdfPath;
        public final String molPath;
        public final String xyzPath;
        public final boolean success;
        public final String message;

        public ExportResult(String name, String sdf, String mol,
                String xyz, boolean ok, String msg) {
            this.name = name;
            this.sdfPath = sdf;
            this.molPath = mol;
            this.xyzPath = xyz;
            this.success = ok;
            this.message = msg;
        }
    }

    // ── SDF export via CDK SDFWriter ──────────────────────────────────────────
    /**
     * Write a molecule to SDF format using CDK's SDFWriter.
     *
     * @param mol molecule with 3D coordinates
     * @param outFile destination .sdf file
     * @return true on success
     */
    public static boolean writeSDFFile(IAtomContainer mol, File outFile) {
        SDFWriter writer = null;
        try {
            writer = new SDFWriter(new FileWriter(outFile));
            writer.write(mol);
            return true;
        } catch (Exception e) {
            System.out.println("    SDF write error: " + e.getMessage());
            return false;
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (Exception ignored) {
                }
            }
        }
    }

    // ── MOL (V2000) export via CDK MDLV2000Writer ────────────────────────────
    /**
     * Write a molecule to MDL MOL V2000 format. MOL is a single-structure SDF —
     * useful for direct import into Avogadro, GaussView, Mercury, and other
     * molecular editors.
     *
     * @param mol molecule with 3D coordinates
     * @param outFile destination .mol file
     * @param title molecule title written to the header line
     * @return true on success
     */
    public static boolean writeMOLFile(IAtomContainer mol, File outFile,
            String title) {
        MDLV2000Writer writer = null;
        try {
            mol.setTitle(title);
            writer = new MDLV2000Writer(new FileWriter(outFile));
            writer.writeMolecule(mol);
            return true;
        } catch (Exception e) {
            System.out.println("    MOL write error: " + e.getMessage());
            return false;
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (Exception ignored) {
                }
            }
        }
    }

    // ── XYZ export (manual — no library needed) ───────────────────────────────
    /**
     * Write a molecule to XYZ format.
     *
     * XYZ format specification: Line 1: number of atoms (integer) Line 2:
     * comment / title (free text) Lines 3+: Element X Y Z (space-separated,
     * Angstrom)
     *
     * This is the primary input format for ORCA and Gaussian. Example ORCA
     * usage: ! B3LYP 6-31G(d) OPT FREQ * xyzfile 0 1 TNT.xyz
     *
     * @param mol molecule with 3D coordinates
     * @param outFile destination .xyz file
     * @param title comment line (line 2 of XYZ file)
     * @return true on success
     */
    public static boolean writeXYZFile(IAtomContainer mol, File outFile,
            String title) {
        FileWriter fw = null;
        try {
            fw = new FileWriter(outFile);
            int n = mol.getAtomCount();

            // Line 1: atom count
            fw.write(n + "\n");

            // Line 2: title/comment
            fw.write(title + " -- MMFF94-optimized geometry, ready for DFT\n");

            // Lines 3+: one atom per line
            for (int i = 0; i < n; i++) {
                IAtom atom = mol.getAtom(i);
                Point3d p = atom.getPoint3d();
                String sym = atom.getSymbol();

                if (p != null) {
                    fw.write(String.format("  %-4s  %14.8f  %14.8f  %14.8f%n",
                            sym, p.x, p.y, p.z));
                } else {
                    fw.write(String.format("  %-4s  %14.8f  %14.8f  %14.8f%n",
                            sym, 0.0, 0.0, 0.0));
                }
            }
            return true;

        } catch (Exception e) {
            System.out.println("    XYZ write error: " + e.getMessage());
            return false;
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (Exception ignored) {
                }
            }
        }
    }

    // ── Combined export: SDF + MOL + XYZ ─────────────────────────────────────
    /**
     * Export a molecule to all three formats in one call.
     *
     * @param name compound name (used for filename and titles)
     * @param mol molecule with 3D coordinates
     * @param outputDir directory to write files into
     * @return ExportResult with paths to all created files
     */
    public static ExportResult exportAll(String name, IAtomContainer mol,
            String outputDir) {
        String safe = name.replace(" ", "_").replace("-", "_");
        File dir = new File(outputDir);
        dir.mkdirs();

        File sdfFile = new File(dir, safe + ".sdf");
        File molFile = new File(dir, safe + ".mol");
        File xyzFile = new File(dir, safe + ".xyz");

        boolean sdfOk = writeSDFFile(mol, sdfFile);
        boolean molOk = writeMOLFile(mol, molFile, name);
        boolean xyzOk = writeXYZFile(mol, xyzFile,
                name + " | SMILES-derived MMFF94 geometry");

        boolean allOk = sdfOk && molOk && xyzOk;
        String status = allOk ? "all formats written"
                : String.format("SDF=%b MOL=%b XYZ=%b", sdfOk, molOk, xyzOk);

        return new ExportResult(name,
                sdfOk ? sdfFile.getPath() : null,
                molOk ? molFile.getPath() : null,
                xyzOk ? xyzFile.getPath() : null,
                allOk, status);
    }

    // ── Print XYZ block to console ────────────────────────────────────────────
    /**
     * Print the full XYZ file content to System.out for inspection. Useful when
     * running inside NetBeans — shows the file content directly in the Output
     * window without needing to open files.
     */
    public static void printXYZToConsole(String name, IAtomContainer mol) {
        int n = mol.getAtomCount();
        System.out.println("  --- XYZ content for " + name + " ---");
        System.out.println("  " + n);
        System.out.println("  " + name + " MMFF94-optimized starting geometry for DFT");
        for (int i = 0; i < n; i++) {
            IAtom atom = mol.getAtom(i);
            Point3d p = atom.getPoint3d();
            if (p != null) {
                System.out.printf("  %-4s  %12.6f  %12.6f  %12.6f%n",
                        atom.getSymbol(), p.x, p.y, p.z);
            }
        }
        System.out.println("  --- end XYZ ---");
    }

    // ── Main ─────────────────────────────────────────────────────────────────
    public static void main(String[] args) {

        System.out.println(rep('=', 62));
        System.out.println("  CHAPTER 2 - STEP 3: Export 3D Structures to File Formats");
        System.out.println("  SDF | MOL | XYZ  -- output to structures/ folder");
        System.out.println("  Based on: Klapötke (2019), Ch. 2 — Theoretical Methods");
        System.out.println(rep('=', 62));
        System.out.println();

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser parser = new SmilesParser(builder);

        // ── Dataset: SMILES for export ────────────────────────────────────────
        Map<String, String> dataset = new LinkedHashMap<String, String>();
        dataset.put("TNT", "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]");
        dataset.put("RDX", "C1N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]");
        dataset.put("HMX", "C1N2CN([N+](=O)[O-])CN([N+](=O)[O-])CN([N+](=O)[O-])C2[N+](=O)[O-]1");
        dataset.put("PETN", "C(CO[N+](=O)[O-])(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]");
        dataset.put("FOX-7", "NC(=C([N+](=O)[O-])[N+](=O)[O-])N");
        dataset.put("NTO", "O=c1[nH]nno1");

        System.out.printf("  %-18s  %-8s  %-8s  %-8s  %s%n",
                "Compound", "SDF", "MOL", "XYZ", "Status");
        System.out.println("  " + rep('-', 65));

        for (Map.Entry<String, String> entry : dataset.entrySet()) {
            String name = entry.getKey();
            String smiles = entry.getValue();

            // Generate 3D structure (reuses Step 1)
            Step1_3DGeneration.Structure3D s
                    = Step1_3DGeneration.generateStructure3D(name, smiles, parser);

            if (!s.success || s.mol == null) {
                System.out.printf("  %-18s  3D generation failed: %s%n",
                        name, s.message);
                continue;
            }

            // Export to all formats
            ExportResult r = exportAll(name, s.mol, "structures");

            System.out.printf("  %-18s  %-8s  %-8s  %-8s  %s%n",
                    name,
                    r.sdfPath != null ? "OK" : "FAIL",
                    r.molPath != null ? "OK" : "FAIL",
                    r.xyzPath != null ? "OK" : "FAIL",
                    r.message);
        }

        System.out.println();
        System.out.println("  Files written to: "
                + new File("structures").getAbsolutePath());

        // ── Print XYZ for nitromethane as a teaching example ─────────────────
        System.out.println();
        System.out.println("  WORKED EXAMPLE — XYZ file content for Nitromethane:");
        System.out.println("  (This is what you paste into an ORCA .inp file)");
        System.out.println();

        Step1_3DGeneration.Structure3D nm
                = Step1_3DGeneration.generateStructure3D(
                        "Nitromethane", "C[N+](=O)[O-]", parser);
        if (nm.success && nm.mol != null) {
            printXYZToConsole("Nitromethane", nm.mol);
        }

        // ── ORCA input snippet ────────────────────────────────────────────────
        System.out.println();
        System.out.println("  ORCA INPUT SNIPPET (using the XYZ file above):");
        System.out.println("  " + rep('-', 50));
        System.out.println("  ! B3LYP 6-31G(d) OPT FREQ TightSCF");
        System.out.println("  %pal nprocs 4 end");
        System.out.println("  %maxcore 2000");
        System.out.println("  * xyzfile 0 1 Nitromethane.xyz");
        System.out.println("  " + rep('-', 50));
        System.out.println("  Reference: Neese et al. (2020). J. Chem. Phys., 152, 224108.");
        System.out.println();

        // ── Format comparison table ───────────────────────────────────────────
        System.out.println("  FORMAT COMPARISON:");
        System.out.println("  " + rep('-', 62));
        System.out.printf("  %-8s  %-10s  %-12s  %-20s%n",
                "Format", "Extension", "Multi-mol?", "Primary Use");
        System.out.println("  " + rep('-', 62));
        System.out.printf("  %-8s  %-10s  %-12s  %-20s%n",
                "SDF", "   .sdf", "   YES", "Database, RDKit, CDK");
        System.out.printf("  %-8s  %-10s  %-12s  %-20s%n",
                "MOL", "   .mol", "   NO ", "Avogadro, GaussView");
        System.out.printf("  %-8s  %-10s  %-12s  %-20s%n",
                "XYZ", "   .xyz", "   NO ", "ORCA, Gaussian, CP2K");
        System.out.println(rep('=', 62));
    }

    static String rep(char c, int n) {
        return new String(new char[n]).replace('\0', c);
    }
}
