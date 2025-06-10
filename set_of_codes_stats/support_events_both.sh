#!/bin/bash

gawk -F"\t" '
{
    split($1, arr, "_");
    gene = arr[1];
    rank = arr[2];
    type = $9;
    col6 = $6;

    if (rank > maximal[gene]) {
        maximal[gene] = rank;
    }

    lines[gene,rank,NR] = $0;
    indices[gene,rank]++;

    if (col6 == "T") {
        if (type == "canonical")
            canonical_exists[gene,rank] = 1;
        else if (type == "alternative")
            alternative_exists[gene,rank] = 1;

        for (i=2; i<=5; i++) {
            tool = i-1;
            if (type == "canonical" && $i == "T")
                canonical_tool_T[gene,rank,tool] = 1;
            if (type == "alternative" && $i == "T")
                alternative_tool_T[gene,rank,tool] = 1;
        }
    }
}
END {
    print "gene,1,2,3,4" > "support_events_both_mouse.csv";

    for (key in canonical_exists) {
        split(key, parts, SUBSEP);
        gene = parts[1];
        rank = parts[2];

        if (alternative_exists[gene,rank]) {
            total_rank_count[gene]++;

            tool_match_count = 0;
            for (tool=1; tool<=4; tool++) {
                if (canonical_tool_T[gene,rank,tool] && alternative_tool_T[gene,rank,tool])
                    tool_match_count++;
            }

            for (line_idx in lines) {
                split(line_idx, p, SUBSEP);
                g = p[1]; r = p[2]; nr = p[3];
                if (g == gene && r == rank) {
                    split(lines[line_idx], fields, "\t");
                    if (fields[6] == "T" && fields[9] == "alternative") {
                        count_1_tool[gene] += (tool_match_count == 1);
                        count_2_tool[gene] += (tool_match_count == 2);
                        count_3_tool[gene] += (tool_match_count == 3);
                        count_4_tool[gene] += (tool_match_count == 4);
                    }
                }
            }
        }
    }

    for (g in total_rank_count) {
        total = total_rank_count[g];
        maxr = maximal[g];
        if (total > 0) {
            printf "%s,%.2f,%.2f,%.2f,%.2f\n", g,
                (count_1_tool[g]/total)*100,
                (count_2_tool[g]/total)*100,
                (count_3_tool[g]/total)*100,
                (count_4_tool[g]/total)*100 >> "support_events_both_mouse.csv";
        }
    }
}
' all_events_mouse.csv
