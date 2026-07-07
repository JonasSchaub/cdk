package org.openscience.cdk.fingerprint;

import java.util.ArrayList;
import java.util.List;

/**
 * Immutable implementation of {@link IntArrayCountFingerprint} storing the occurrence
 * count of each Biosynfoni substructure.
 * <p>
 * The fingerprint is created from the filtered SMARTS matches generated during
 * fingerprint calculation. Each fingerprint position corresponds to one
 * Biosynfoni substructure key, and the stored value represents the number of matches for that key.
 * <p>
 * This implementation provides read-only access to the fingerprint counts.
 * Methods intended for mutable fingerprints are currently not supported.
 */
final class BiosynfoniCountFingerprint extends IntArrayCountFingerprint {
    public BiosynfoniCountFingerprint(List<List<int[]>> filteredMatches) {
        int[] counts = new int[filteredMatches.size()];

        for (int i = 0; i < filteredMatches.size(); i++) {
            counts[i] = filteredMatches.get(i).size();
        }

        List<Integer> hashes = new ArrayList<>();
        List<Integer> values = new ArrayList<>();

        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > 0) {
                hashes.add(i);
                values.add(counts[i]);
            }
        }

        this.hitHashes = new int[hashes.size()];
        this.numOfHits = new int[values.size()];

        for (int i = 0; i < hashes.size(); i++) {
            this.hitHashes[i] = hashes.get(i);
            this.numOfHits[i] = values.get(i);
        }
    }
}
