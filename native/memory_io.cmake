# Adapt generated copies, leaving the pinned vendor checkout unchanged. Limit
# replacements to the writer functions; fail configuration if their anchors
# change rather than silently falling back to disk I/O after a vendor upgrade.
function(pf_replace_required variable before after)
    string(FIND "${${variable}}" "${before}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "IQ-TREE memory adapter anchor missing: ${before}")
    endif()
    string(REPLACE "${before}" "${after}" modified "${${variable}}")
    set(${variable} "${modified}" PARENT_SCOPE)
endfunction()

function(pf_memory_source target relative begin end)
    set(original "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/iqtree2/${relative}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${original}")
    file(READ "${original}" source)
    string(FIND "${source}" "${begin}" first)
    string(FIND "${source}" "${end}" last)
    if(first EQUAL -1 OR last EQUAL -1 OR NOT last GREATER first)
        message(FATAL_ERROR "IQ-TREE writer boundaries changed: ${relative}")
    endif()
    math(EXPR length "${last} - ${first}")
    string(SUBSTRING "${source}" 0 ${first} prefix)
    string(SUBSTRING "${source}" ${first} ${length} writer)
    string(SUBSTRING "${source}" ${last} -1 suffix)
    pf_replace_required(writer "new ofstream(" "new pf::MemoryOutputFile(")
    pf_replace_required(writer "((ofstream*)" "((pf::MemoryOutputFile*)")
    set(source "${prefix}${writer}${suffix}")
    # Numeric mode bypasses FASTA and indel numeric-text serialization.
    if(target STREQUAL "simulator")
        pf_replace_required(source
            "void AliSimulator::writeInternalStatesIndels(Node* node, ostream &out)\n{"
            "void AliSimulator::writeInternalStatesIndels(Node* node, ostream &out)\n{\n    if (pf::raw_enabled()) {\n        pf::raw_state().originals.push_back({node->name, std::move(node->sequence->sequence_chunks[0])});\n        map_seqname_node[node->name] = node;\n        return;\n    }")
        pf_replace_required(source
            "    // we can only write and delete sequence chunk in normal simulations: without Indels, FunDi, ASC, etc"
            "    if (pf::raw_enabled() && store_seq_at_cache) {\n        if ((*it)->node->isLeaf()) pf::capture_raw((*it)->node->name, node_seq_chunk, segment_length, state_mapping);\n        if (node->isLeaf() && node->name != ROOT_NAME) pf::capture_raw(node->name, dad_seq_chunk, segment_length, state_mapping);\n        return;\n    }\n    // we can only write and delete sequence chunk in normal simulations: without Indels, FunDi, ASC, etc")
    else()
        pf_replace_required(source
            "    // read tmp_data line by line"
            "    if (pf::raw_enabled()) { pf::capture_indel_rows(alisimulator, sequence_length, state_mapping); return; }\n    // read tmp_data line by line")
    endif()
    if(target STREQUAL "main")
        pf_replace_required(source
            "    // if write_sequences_from_tmp_data and this node is a leaf -> skip this node as its sequence was already written to the output file"
            "    if (pf::raw_enabled()) { pf::capture_retained_rows(node, dad, write_sequences_from_tmp_data, sequence_length, state_mapping); return; }\n    // if write_sequences_from_tmp_data and this node is a leaf -> skip this node as its sequence was already written to the output file")
        pf_replace_required(source "igzstream in;" "pf::MemoryInputFile in;")
        pf_replace_required(source "An alignment written to " "An alignment returned in memory for ")
        pf_replace_required(source "Unaligned sequences written to " "Gap-free sequences returned in memory for ")
        pf_replace_required(source
            "remove((super_alisimulator->params->alisim_output_filename + \"_\" + super_alisimulator->params->tmp_data_filename + \"_\" + convertIntToString(MPIHelper::getInstance().getProcessID())).c_str());"
            "pf::memory_files().erase(super_alisimulator->params->alisim_output_filename + \"_\" + super_alisimulator->params->tmp_data_filename + \"_\" + convertIntToString(MPIHelper::getInstance().getProcessID()));")
    endif()

    get_filename_component(filename "${original}" NAME)
    get_filename_component(source_directory "${original}" DIRECTORY)
    set(generated "${CMAKE_CURRENT_BINARY_DIR}/memory_${filename}")
    set(content "#include \"${CMAKE_CURRENT_SOURCE_DIR}/native/memory_io.hpp\"\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/native/numeric_output.hpp\"\n${source}")
    # configure_file avoids touching the copy when its contents are unchanged.
    set(template "${CMAKE_CURRENT_BINARY_DIR}/memory_${filename}.in")
    file(WRITE "${template}" "@content@")
    configure_file("${template}" "${generated}" @ONLY)
    get_target_property(sources ${target} SOURCES)
    list(REMOVE_ITEM sources "${filename}")
    list(APPEND sources "${generated}")
    set_property(TARGET ${target} PROPERTY SOURCES "${sources}")
    target_include_directories(${target} PRIVATE "${source_directory}")
endfunction()

pf_memory_source(main "main/alisim.cpp"
    "void writeSequencesToFile(" "void mergeAndWriteSequencesToFiles(")
pf_memory_source(simulator "simulator/alisimulator.cpp"
    "void AliSimulator::openOutputStream(" "void AliSimulator::simulateSeqs(int")
