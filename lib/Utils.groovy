//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    // Check if database filepaths exists
    public static Boolean fileExists(filename) {
        return new File(filename).exists()
    }

     //
    // Load metadata from TSV file and apply filters
    // Returns a map of sample_ID -> metadata_map
    //
    public static Map loadMetadataWithFilters(metadata_file, filter_map, log) {
        def metadata_map = [:]
        def metadata_lines = metadata_file.readLines()
        
        if (metadata_lines.size() == 0) {
            log.warn "Metadata file is empty"
            return metadata_map
        }
        
        // Parse header
        def header = metadata_lines[0].split('\t')
        def col_indices = [:]
        header.eachWithIndex { col, idx -> col_indices[col.trim()] = idx }
        
        // Required columns
        def required_cols = ['sample_ID']
        required_cols.each { col ->
            if (!col_indices.containsKey(col)) {
                throw new Exception("Metadata file missing required column: ${col}")
            }
        }
        
        // Validate that filter columns exist in metadata
        filter_map.each { filter_col, filter_val ->
            if (!col_indices.containsKey(filter_col)) {
                log.warn "Filter column '${filter_col}' not found in metadata - ignoring this filter"
            }
        }
        
        // Process data rows
        metadata_lines[1..-1].each { line ->
            if (line.trim() == '' || line.startsWith('#')) {
                return
            }
            
            def fields = line.split('\t')
            def sample_id = fields[col_indices['sample_ID']].trim()
            
            // Apply filters - all must match (AND logic)
            def passes_filters = true
            
            filter_map.each { filter_col, filter_val ->
                if (col_indices.containsKey(filter_col)) {
                    def col_idx = col_indices[filter_col]
                    if (col_idx < fields.size()) {
                        def actual_val = fields[col_idx].trim()
                        if (actual_val != filter_val) {
                            passes_filters = false
                        }
                    } else {
                        passes_filters = false
                    }
                }
            }
            
            if (passes_filters) {
                // Store all metadata for this sample
                def sample_metadata = [:]
                col_indices.each { col_name, idx ->
                    if (idx < fields.size()) {
                        sample_metadata[col_name] = fields[idx].trim()
                    }
                }
                metadata_map[sample_id] = sample_metadata
            }
        }
        
        return metadata_map
    }
}