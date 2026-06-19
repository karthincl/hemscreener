//Practical Chemoinformatics,2014, M.Karthikeyan, Renu Vyas;

package chemoinformatics;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// CDK imports for molecule validation
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class Step3_NameToSMILES {

    private static final String PUBCHEM_BASE
            = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/";

    private static final String PROPERTIES
            = "IsomericSMILES,CanonicalSMILES,MolecularFormula,MolecularWeight,InChI,InChIKey";

    public static class CompoundInfo {

        public final String canonicalSMILES;
        public final String molecularFormula;
        public final String molecularWeight;
        public final String inchi;
        public final String inchiKey;
        public final int cid;

        public CompoundInfo(String canonical, String formula, String mw, String inchi, String inchiKey, int cid) {
            this.canonicalSMILES = canonical;
            this.molecularFormula = formula;
            this.molecularWeight = mw;
            this.inchi = inchi;
            this.inchiKey = inchiKey;
            this.cid = cid;
        }
    }

    public static CompoundInfo parsePubChemJson(String body) {
        // Updated keys matching your raw JSON layout string properties precisely
        String smiles = extractJsonString(body, "SMILES");
        String formula = extractJsonString(body, "MolecularFormula");
        String mw = extractJsonString(body, "MolecularWeight");
        String inchi = extractJsonString(body, "InChI");
        String inchiKey = extractJsonString(body, "InChIKey");
        int cid = extractCID(body);

        if (smiles == null || smiles.isEmpty()) {
            return null;
        }

        return new CompoundInfo(smiles, formula, mw, inchi, inchiKey, cid);
    }

    private static String extractJsonString(String json, String key) {

        String pattern = "\"" + key + "\":";
        int start = json.indexOf(pattern);
        if (start < 0) {
            return null;
        }

        start += pattern.length();

        while (start < json.length() && (json.charAt(start) == ' ' || json.charAt(start) == '"')) {
            start++;
        }

        int end = json.indexOf("\"", start);
        if (end < 0) {
            return null;
        }
        return json.substring(start, end);
    }

    private static int extractCID(String json) {
        String pattern = "\"CID\":";
        int idx = json.indexOf(pattern);
        if (idx < 0) {
            return -1;
        }

        idx += pattern.length();
        while (idx < json.length() && (json.charAt(idx) == ' ' || json.charAt(idx) == ':')) {
            idx++;
        }

        int end = idx;
        while (end < json.length() && Character.isDigit(json.charAt(end))) {
            end++;
        }

        try {
            return Integer.parseInt(json.substring(idx, end));
        } catch (NumberFormatException e) {
            return -1;
        }
    }

    static String moldata = "";
    static String failed_mol = "";

    static IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
    static SmilesParser parser = new SmilesParser(builder);

    public static CompoundInfo nameToSMILES(String compoundName) {

        HttpURLConnection connection = null;
        try {

            String encoded = URLEncoder.encode(compoundName, StandardCharsets.UTF_8.name());
            String urlStr = PUBCHEM_BASE + encoded + "/property/" + PROPERTIES + "/JSON";
            URL url = new URL(urlStr);
            connection = (HttpURLConnection) url.openConnection();
            connection.setConnectTimeout(15000); // 15 seconds
            connection.setReadTimeout(20000);    // 20 seconds
            int responseCode = connection.getResponseCode();

            if (responseCode == 404) {
                //System.out.printf("  Not found in PubChem: %s%n", compoundName);
                failed_mol += compoundName + "\n";
                return null;
            }
            if (responseCode != 200) {
                System.out.printf("  HTTP %d for: %s%n", responseCode, compoundName);
                return null;
            }

            StringBuilder responseBody = new StringBuilder();
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(connection.getInputStream(), StandardCharsets.UTF_8))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    responseBody.append(line);
                }
            }
            // System.out.println(responseBody);
            String body = responseBody.toString();
            CompoundInfo info = parsePubChemJson(body);
            if (info != null) {

                moldata += compoundName + "\t";
                moldata += info.cid + "\t";
                moldata += info.molecularFormula + "\t";
                moldata += info.molecularWeight + "\t";
                moldata += info.canonicalSMILES + "\t";
                moldata += info.inchi + "\t";
                moldata += info.inchiKey + "\t";

                try {
                    IAtomContainer mol = parser.parseSmiles(info.canonicalSMILES);
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
                    moldata += mol.getAtomCount() + "\t";
                    moldata += mol.getBondCount() + "\n";

                } catch (Exception e) {

                }
            } else {
                System.out.println("Extraction failed! SMILES data field not located.");
            }
            return new CompoundInfo(info.canonicalSMILES, info.molecularFormula, info.molecularWeight, info.inchi, info.inchiKey, info.cid);
        } catch (Exception e) {
            return null;
        } finally {
            if (connection != null) {
                connection.disconnect();
            }
        }
    }

    private static String extractField(String source, String regex) {
        Pattern pattern = Pattern.compile(regex);
        Matcher matcher = pattern.matcher(source);
        if (matcher.find()) {
            return matcher.group(1);
        }
        return "?";
    }

    private static String extractJsonValue(String json, String key) {
        String pattern = "\"" + key + "\":";
        int start = json.indexOf(pattern);
        if (start < 0) {
            return null;
        }
        start += pattern.length();
        if (start < json.length() && json.charAt(start) == '"') {
            start++;
            int end = json.indexOf("\"", start);
            return end < 0 ? null : json.substring(start, end);
        }
        int end = start;
        while (end < json.length()
                && (Character.isDigit(json.charAt(end)) || json.charAt(end) == '.')) {
            end++;
        }
        return json.substring(start, end);
    }

    private static String generateRepeatedChar(char ch, int count) {
        char[] array = new char[count];
        java.util.Arrays.fill(array, ch);
        return new String(array);
    }

    public static void main(String[] args) {

        System.out.println("============================================================");
        System.out.println("  Name to SMILES via PubChem REST API (Java 1.8)");
        System.out.println("  Requires: active internet connection");
        System.out.println("============================================================\n");
        disableSSLVerification();

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
            "FOX-7",};

        for (String name : names) {
            try {
                Thread.sleep(250);
            } catch (InterruptedException ignored) {
            }
            CompoundInfo info = nameToSMILES(name);
        }

        System.out.print("Compound Name\tCID\tFormula\tMol Wt\tSMILES\tInChI\tInChIKey\tAtom count\tBond Count\n" + moldata);

        System.out.println("\n\n======Not found in PubChem:=====\n" + failed_mol);

        System.out.println("\nCitation: \n1) Practical Chemoinformatics,2014, M.Karthikeyan, Renu Vyas; \n2) Kim et al. (2023). PubChem 2023 update.");
        System.out.println("  Nucleic Acids Res., 51, D1373-D1380.");
        System.out.println("======= ");
    }

    // Add this method inside your class
    private static void disableSSLVerification() {
        try {
            javax.net.ssl.TrustManager[] trustAllCerts = new javax.net.ssl.TrustManager[]{
                new javax.net.ssl.X509TrustManager() {
                    public java.security.cert.X509Certificate[] getAcceptedIssuers() {
                        return null;
                    }

                    public void checkClientTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                    }

                    public void checkServerTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                    }
                }
            };

            javax.net.ssl.SSLContext sc = javax.net.ssl.SSLContext.getInstance("SSL");
            sc.init(null, trustAllCerts, new java.security.SecureRandom());
            javax.net.ssl.HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());

            // Also bypass hostname verification
            javax.net.ssl.HttpsURLConnection.setDefaultHostnameVerifier((hostname, session) -> true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
