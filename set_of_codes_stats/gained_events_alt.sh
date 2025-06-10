gawk -F"\t" '
{
    split($1, arr, "_");
    gene = arr[1];
    rank = arr[2];
    type = $9;
    col6 = $6;

    # Save info per line
    lines[gene,rank,NR] = $0;
    line_ids[gene,rank][NR] = 1;

    # Save canonical T and alternative F
    if (type == "canonical" && col6 == "T") {
        canonical_t[gene,rank] = 1;
    } else if (type == "alternative" && col6 == "F") {
        alternative_f[gene,rank] = 1;

        # Count tools only for alternative with col6 == F
        tool_count = 0;
        for (i = 2; i <= 5; i++)
            if ($i == "T") tool_count++;
        tools[gene,rank,NR] = tool_count;
    }
}
END {
    print "gene,1,2,3,4" > "gained_events_alt_mouse.csv";

    for (key in canonical_t) {
        if (alternative_f[key]) {
            split(key, parts, SUBSEP);
            gene = parts[1];
            rank = parts[2];

            total_rank_count[gene]++;

            # Go through saved line IDs for this gene+rank
            for (nr in line_ids[gene,rank]) {
                split(lines[gene,rank,nr], fields, "\t");
                if (fields[9] == "alternative" && fields[6] == "F") {
                    count = tools[gene,rank,nr];
                    if (count == 1) count_1[gene]++;
                    else if (count == 2) count_2[gene]++;
                    else if (count == 3) count_3[gene]++;
                    else if (count == 4) count_4[gene]++;
                }
            }
        }
    }

    for (g in total_rank_count) {
        total = total_rank_count[g];
        if (total > 0) {
            printf "%s,%.2f,%.2f,%.2f,%.2f\n", g,
                (count_1[g]/total)*100,
                (count_2[g]/total)*100,
                (count_3[g]/total)*100,
                (count_4[g]/total)*100 >> "gained_events_alt_mouse.csv";
        }
    }
}
' all_events_mouse.csv
