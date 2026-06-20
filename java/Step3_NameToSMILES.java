package chemoinformatics;

// ============================================================
//  Step 3 — Name-to-SMILES via PubChem REST API
//  Chapter 1: Introduction and Environment Setup
//  Based on: Chemistry of High-Energy Materials, 4th ed.
//             Klapötke, T. M. (De Gruyter, 2019)
// ------------------------------------------------------------
//  Equivalent to Python:  step4_name_to_smiles.py  (uses pubchempy)
//  Java uses the PubChem REST API directly via java.net.http.HttpClient
//  (built into Java 11+ — no external library needed)
//
//  Run in NetBeans: Right-click -> Run File (Shift+F6)
//  Note: requires internet connection to reach pubchem.ncbi.nlm.nih.gov
// ============================================================

import java.net.URI;
import java.net.URLEncoder;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.LinkedHashMap;
import java.util.Map;

// CDK imports for molecule validation
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.interfaces.IMolecularFormula;

/**
 * Step 3 — Convert Compound Names to SMILES Using the PubChem REST API.
 *
 * In Python, pubchempy.get_compounds() handles this automatically.
 * In Java, we call the PubChem REST API directly using Java 11's
 * built-in HttpClient (no external dependencies required).
 *
 * PubChem REST API endpoint used:
 *   https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{NAME}/property/
 *       IsomericSMILES,CanonicalSMILES,MolecularFormula,
 *       MolecularWeight,InChI,InChIKey/JSON
 *
 * The JSON response is parsed with simple string operations to avoid
 * adding a JSON library dependency. In production code, use Gson or Jackson.
 *
 * Citation:
 *   Kim, S. et al. (2023). PubChem 2023 update.
 *   Nucleic Acids Research, 51(D1), D1373-D1380.
 *   https://doi.org/10.1093/nar/gkac956
 */
public class Step3_NameToSMILES {

    // Base URL for PubChem PUG REST API
    private static final String PUBCHEM_BASE =
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/";

    // Properties to retrieve in a single API call
    private static final String PROPERTIES =
        "IsomericSMILES,CanonicalSMILES,MolecularFormula,MolecularWeight,InChI,InChIKey";

    // HTTP client — reused across all requests (Java 11+)
    private static final HttpClient HTTP_CLIENT = HttpClient.newBuilder()
        .connectTimeout(Duration.ofSeconds(15))
        .build();

    /**
     * Holds the result of a PubChem lookup.
     * Mirrors the dict returned by name_to_smiles() in the Python tutorial.
     */
    public static class CompoundInfo {
        public final String name;
        public final String canonicalSMILES;
        public final String isomericSMILES;
        public final String molecularFormula;
        public final String molecularWeight;
        public final String inchi;
        public final String inchiKey;
        public final int    cid;

        public CompoundInfo(String name, String canonical, String isomeric,
                            String formula, String mw, String inchi,
                            String inchiKey, int cid) {
            this.name             = name;
            this.canonicalSMILES  = canonical;
            this.isomericSMILES   = isomeric;
            this.molecularFormula = formula;
            this.molecularWeight  = mw;
            this.inchi            = inchi;
            this.inchiKey         = inchiKey;
            this.cid              = cid;
        }
    }

    // ── Core lookup method ────────────────────────────────────────────────────
    /**
     * Query PubChem for a compound by its common name.
     *
     * @param compoundName  e.g. "TNT", "RDX", "ammonium nitrate"
     * @return CompoundInfo  or  null if not found / network error
     */
    public static CompoundInfo nameToSMILES(String compoundName) {
        try {
            // URL-encode the compound name (handles spaces and special chars)
            String encoded = URLEncoder.encode(compoundName, StandardCharsets.UTF_8);
            String url     = PUBCHEM_BASE + encoded + "/property/" + PROPERTIES + "/JSON";

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .timeout(Duration.ofSeconds(20))
                .header("Accept", "application/json")
                .GET()
                .build();

            HttpResponse<String> response =
                HTTP_CLIENT.send(request, HttpResponse.BodyHandlers.ofString());

            if (response.statusCode() == 404) {
                System.out.printf("  Not found in PubChem: %s%n", compoundName);
                return null;
            }
            if (response.statusCode() != 200) {
                System.out.printf("  HTTP %d for: %s%n", response.statusCode(), compoundName);
                return null;
            }

            // ── Minimal JSON parsing using String operations ──────────────────
            // Avoids requiring Gson/Jackson. In production code use a JSON library.
            String body = response.body();

            String canonical   = extractJsonString(body, "CanonicalSMILES");
            String isomeric    = extractJsonString(body, "IsomericSMILES");
            String formula     = extractJsonString(body, "MolecularFormula");
            String mw          = extractJsonValue(body,  "MolecularWeight");
            String inchi       = extractJsonString(body, "InChI");
            String inchiKey    = extractJsonString(body, "InChIKey");
            int    cid         = extractCID(body);

            if (canonical == null || canonical.isEmpty()) {
                System.out.printf("  No SMILES returned for: %s%n", compoundName);
                return null;
            }

            return new CompoundInfo(compoundName, canonical, isomeric,
                                    formula, mw, inchi, inchiKey, cid);

        } catch (java.io.IOException | InterruptedException e) {
            System.out.printf("  Network error for '%s': %s%n",
                              compoundName, e.getMessage());
            return null;
        }
    }

    // ── JSON extraction helpers ───────────────────────────────────────────────

    /** Extract a quoted string value for a given JSON key. */
    private static String extractJsonString(String json, String key) {
        String pattern = "\"" + key + "\":\"";
        int start = json.indexOf(pattern);
        if (start < 0) return null;
        start += pattern.length();
        int end = json.indexOf("\"", start);
        if (end < 0) return null;
        return json.substring(start, end);
    }

    /** Extract an unquoted numeric value for a given JSON key. */
    private static String extractJsonValue(String json, String key) {
        String pattern = "\"" + key + "\":";
        int start = json.indexOf(pattern);
        if (start < 0) return null;
        start += pattern.length();
        // Skip optional quote
        if (start < json.length() && json.charAt(start) == '"') {
            start++;
            int end = json.indexOf("\"", start);
            return end < 0 ? null : json.substring(start, end);
        }
        int end = start;
        while (end < json.length() &&
               (Character.isDigit(json.charAt(end)) || json.charAt(end) == '.')) {
            end++;
        }
        return json.substring(start, end);
    }

    /** Extract the first CID (compound identifier) from the JSON. */
    private static int extractCID(String json) {
        String pattern = "\"CID\":";
        int idx = json.indexOf(pattern);
        if (idx < 0) return -1;
        idx += pattern.length();
        int end = idx;
        while (end < json.length() && Character.isDigit(json.charAt(end))) end++;
        try {
            return Integer.parseInt(json.substring(idx, end));
        } catch (NumberFormatException e) {
            return -1;
        }
    }

    // ── Main: run lookup for energetic compounds ──────────────────────────────

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  STEP 3 — Name to SMILES via PubChem REST API");
        System.out.println("  Requires: active internet connection");
        System.out.println("============================================================\n");

        // Energetic compounds from Klapötke (2019), Chapter 1
        String[] names = {
            "TNT",
            "RDX",
            "HMX",
            "PETN",
            "nitromethane",
            "ammonium nitrate",
            "picric acid",
            "tetryl",
            "TATB",
            "FOX-7",
        };

        System.out.printf("  %-22s  %6s  %-12s  %8s  %s%n",
                          "Name", "CID", "Formula", "MW", "Canonical SMILES (truncated)");
        System.out.println("  " + "-".repeat(90));

        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser       parser  = new SmilesParser(builder);

        for (String name : names) {
            // Add a small delay to be a polite API client (PubChem rate limit ~5 req/s)
            try { Thread.sleep(250); } catch (InterruptedException ignored) {}

            CompoundInfo info = nameToSMILES(name);

            if (info == null) {
                System.out.printf("  %-22s  %6s  %-12s  %8s  %s%n",
                                  name, "?", "?", "?", "(not found)");
                continue;
            }

            // Truncate SMILES for display
            String smiDisplay = info.canonicalSMILES.length() > 45
                ? info.canonicalSMILES.substring(0, 45) + "..."
                : info.canonicalSMILES;

            System.out.printf("  %-22s  %6d  %-12s  %8s  %s%n",
                              info.name, info.cid, info.molecularFormula,
                              info.molecularWeight, smiDisplay);

            // ── Validate: round-trip through CDK ──────────────────────────────
            // Parse the SMILES returned by PubChem using CDK to confirm validity
            try {
                IAtomContainer mol = parser.parseSmiles(info.canonicalSMILES);
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                // If we reach here, CDK accepted the SMILES as valid
            } catch (Exception e) {
                System.out.printf("    WARNING: CDK could not parse returned SMILES: %s%n",
                                  e.getMessage());
            }
        }

        System.out.println("\n  TIP: Store results in a Map<String,CompoundInfo> for");
        System.out.println("  subsequent descriptor calculations (Step 5 onwards).");
        System.out.println("\n  Citation: Kim et al. (2023). PubChem 2023 update.");
        System.out.println("            Nucleic Acids Res., 51, D1373-D1380.");
        System.out.println("============================================================");
    }
}
