import os
import matplotlib.pyplot as plt
import pathlib as path
import numpy as np
import datetime


# Get the paths of all elements from a directory
def get_file_paths(directory):
    return [f for f in path.Path(directory).iterdir() if f.is_file()]

# a & b are list of tuples (txt, float) ...
# sum the second component, keep the first one taken from a (must be in same order)
def sum_tuple(a, b):
    a = list(a)
    b = list(b)
    c = []
    for i in range(len(a)):
        c.append((a[i][0], a[i][1]+b[i][1]))
    c = tuple(c)
    return c


# Read all the files in a directory and create a dictionary result
def read_all(directories):
    results = {}
    check_results = {}
    bad_runs = {}
    for dir in directories:
        for file_path in get_file_paths(dir):
            fname = file_path.name
            c = fname.split("-")  # "components"
            if len(c) == 5:
                c[2] = c[2][:-4]  # remove ".txt"
                slen = int(c[3])
                wratio = float(c[4])
                run_name = "-".join(c[1:])  # All to check that all exec yield the same result for a given run
                execdic = results.setdefault(c[0], {})
                datadic = execdic.setdefault(c[1], {})
                datadic.setdefault("nb", int(0))
                datadic.setdefault("computation_ratios", (("lb_Kim", 0), ("lb_Keogh1", 0), ("lb_Keogh2", 0), ("DTW", 0)))
                querydic = datadic.setdefault(c[2], {})
                lendic = querydic.setdefault(slen, {})
                ### Read file
                with open(file_path) as f:
                    l = f.readlines()
                    ### Checking
                    loc = (l[1].split(':')[1]).strip()
                    dist = (l[2].split(':')[1]).strip()
                    v = ",".join([loc, dist])
                    res = check_results.setdefault(run_name, v)
                    if res != v: bad_runs[run_name] = (c[0], c[1], c[2], slen, wratio)
                    ### Time
                    exectime = float((l[4].split(':')[1]).strip().split()[0])
                    lendic[wratio] = (exectime, loc, dist)
                    ### Computation Ratios
                    lbkim = float((l[6].split(':')[1]).strip()[:-1])
                    lbkeogh1 = float((l[7].split(':')[1]).strip()[:-1])
                    lbkeogh2 = float((l[8].split(':')[1]).strip()[:-1])
                    dtw = float((l[9].split(':')[1]).strip()[:-1])
                    datadic["computation_ratios"] = sum_tuple(datadic["computation_ratios"], (("lb_Kim", lbkim), ("lb_Keogh1", lbkeogh1), ("lb_Keogh2", lbkeogh2), ("DTW", dtw)))
                    datadic["nb"] += 1
    # Check results
    for k, v in bad_runs.items():
        print("Inconsistent results for ", k, ":")
        for e in results.keys():
            print("  ", e, results[e][v[1]][v[2]][v[3]][v[4]])
        print("")
    #
    return results


if __name__ == "__main__":

    # Create output directory
    if not os.path.exists('figures'):
        os.makedirs('figures')

    ### Extra dataset info (length)
    datasetInfo = {}
    datasetInfo["FOG_LEG"]   = 1724584
    datasetInfo["SoccerPos"] = 1998606
    datasetInfo["PAMAP2"]    = 3657033
    datasetInfo["MITDB"]     = 27950000
    datasetInfo["REFIT"]     = 78596631
    datasetInfo["PPG"]       = 333570000

    data = read_all(["results/UCR", "results/UCR_USP", "results/UCR_MON", "results/UCR_MON_nolb"])

    # Recover the executables, datasets, queries, lengths and ratios
    # execs = list(data.keys())
    # execs.sort()
    # Force the order:
    execs = ['UCR', 'UCR_USP', 'UCR_MON', 'UCR_MON_nolb']

    # print(execs)
    # Get the list of dataset
    # datasets = list(data[execs[0]].keys())
    # datasets.sort()
    # Force the order: get the datasets sorted by the length of the reference series
    datasets = ['FOG_LEG', 'SoccerPos', 'PAMAP2', 'MITDB', 'REFIT', 'PPG']

    # Name matching original paper
    dataset_names = {'FOG_LEG':'FoG', 'SoccerPos':'Soccer', 'PAMAP2':'PAMAP2', 'MITDB':'ECG', 'REFIT':'REFIT', 'PPG':'PPG'}

    # Queries
    queries = list(q for q in data[execs[0]][datasets[0]].keys() if q.startswith("Query"))
    queries.sort()

    # Lengths
    lengths = list(data[execs[0]][datasets[0]][queries[0]].keys())
    lengths.sort()

    # Window ratios
    ratios = list(data[execs[0]][datasets[0]][queries[0]][lengths[0]].keys())
    ratios.sort()



    ### Extra text per dataset
    dextras = {}
    for d in datasets:
        # Length of the data
        cdtw = "Computed by DTW"
        dl = "Ref length"
        dlnum = f"{datasetInfo[d]:,}"
        dlnum_l = len(dlnum)
        m = len(cdtw)
        # Create the lines of text
        l = ["Proportion pruned by:"]
        # Computation ratios
        for lb, r in data["UCR"][d]["computation_ratios"]:
            s = ""
            if lb == "DTW":
                s = cdtw
            else:
                s = "   " + lb
            avg = r/data["UCR"][d]["nb"]
            num = f"{avg:5.2f}"+"%"
            diff = m - len(s)+2
            l.append(s+(": ".ljust(diff))+num)
        prevlen = len(l[-1])
        l.append(dl+(": ").ljust( prevlen-dlnum_l-len(dl)  )+dlnum)

        # Create the final text per dataset
        dextras[d] = "\n".join(l)
    #print(dextras)

    # Figure base size
    w = 4

    # Same as Fig 10, section 5.2
    # in Speeding up similarity search under dynamic time warping by pruning unpromising alignments
    # * Datasets sorted by length of reference series
    # * All exec on one figure
    # * y axis = average time in seconds
    # * x axis = length of the queries
    # * plotted value = average over queries and window ratios

    WIDTH = 3
    HEIGHT = 2
    FIGWIDTH = WIDTH * w + 1
    FIGHEIGHT = HEIGHT * w + 1
    fig = plt.figure(figsize=(FIGWIDTH, FIGHEIGHT))  # X*Y
    axs_ = fig.subplots(HEIGHT, WIDTH) # nbrow * nbcol
    axstab = []
    for i in range(0, HEIGHT):
        for j in range(0, WIDTH):
            axstab.append(axs_[i, j])

    i = int(0)
    result1024 = {}
    for d in datasets:
        axs = axstab[i]
        i += 1
        axs.set_xlabel('Queries length')
        axs.set_ylabel('Time (s)')
        axs.set_title(dataset_names[d])
        for e in execs:
            result_exec = result1024.setdefault(d, {})
            points = {}
            for q in queries:
                for l in lengths:
                    for r in ratios:
                        times = points.setdefault(l, [])
                        times.append(data[e][d][q][l][r][0])
            y = []
            for (k, v) in points.items():
                y.append(sum(v) / len(v))

            # Save result for l=1024
            result_exec[e] = y[3]

            x = np.arange(len(lengths))  # plot on [0, 1, 2, 3] : evenly spaced
            axs.plot(x, y, label=e)
            axs.set_xticks(x) # Limit ticks to [0, 1, 2, 3]
            axs.set_xticklabels(lengths) # Change the labels to actual lenghts
            axs.ticklabel_format(axis='y', style='scientific', scilimits=[0,4])
            axs.set_ylim(ymin=0) # start plot at 0
            axs.text(0.03, 0.97, dextras[d],
                     ha='left', va='top', transform=axs.transAxes,
                     bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.35'),
                     fontdict={'family': 'monospace'}
                     )
    # End of for d in datasets
    handles, labels = axstab[0].get_legend_handles_labels()
    fig.set_tight_layout({"rect": (0, 0.1, 1, 1.025)})
    leg = fig.legend(handles, labels, ncol=4, loc="upper center", bbox_to_anchor=[0.5, 0.1])
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    #fig.suptitle("For each dataset, average running times per queries' length of all (query, window ratio) combinations", x=0.5, y=0.06)
    #fig.show()
    fig.savefig("figures/time_per_qlength.pdf", format="pdf", bbox_inches = 'tight', pad_inches = 0)



    # Same as figure 11 section 5.2
    # in Speeding up similarity search under dynamic time warping by pruning unpromising alignments
    # * Datasets sorted by length of reference series
    # * All exec on one figure
    # * y axis = average time in seconds
    # * x axis = windows ratios
    # * plotted value = average over queries and queries length

    fig = plt.figure(figsize=(FIGWIDTH, FIGHEIGHT))  # X*Y
    axs_ = fig.subplots(HEIGHT, WIDTH) # nbrow * nbcol
    axstab = []
    for i in range(0, HEIGHT):
        for j in range(0, WIDTH):
            axstab.append(axs_[i, j])

    i = int(0)
    for d in datasets:
        axs = axstab[i]
        i += 1
        axs.set_xticks(ratios)
        axs.set_xlabel('Window ratios')
        axs.set_ylabel('Time (s)')
        axs.set_title(dataset_names[d])
        for e in execs:
            points = {}
            for q in queries:
                for l in lengths:
                    for r in ratios:
                        times = points.setdefault(r, [])
                        times.append(data[e][d][q][l][r][0])
            x = []
            y = []
            for (k, v) in points.items():
                x.append(k)
                y.append(sum(v) / len(v))
            #
            axs.plot(x, y, label=e)
            axs.ticklabel_format(axis='y', style='scientific', scilimits=[0, 4])
            axs.set_ylim(bottom=0)
            axs.text(0.03, 0.97, dextras[d],
                     ha='left', va='top', transform=axs.transAxes,
                     bbox=dict(facecolor='none', edgecolor='black', boxstyle='square,pad=0.35'),
                     fontdict={'family': 'monospace'}
                     )
    # End of for d in datasets
    handles, labels = axstab[0].get_legend_handles_labels()
    fig.set_tight_layout({"rect": (0, 0.1, 1, 1.025)})
    leg = fig.legend(handles, labels, ncol=4, loc="upper center", bbox_to_anchor=[0.5, 0.1])
    for line in leg.get_lines():
        line.set_linewidth(4.0)
     #fig.suptitle("For each dataset, average running times per window ratio of all (query, query's length) combinations", x=0.5, y=0.06)
    #fig.show()
    fig.savefig("figures/time_per_wratio.pdf", format="pdf", bbox_inches = 'tight', pad_inches = 0)



    # Check results of length 1024
    print("--------------------------------------------")
    print("Average timings and speedup for length 1024:")
    for k, v in result1024.items():
        ucr = v["UCR"]
        usp = v["UCR_USP"]
        mon = v["UCR_MON"]
        print(k, f"UCR: {ucr:.2f}s,   UCR USP: {usp:.2f}s,   UCR MON: {mon:.2f}s")
        print(k, f"UCR speedup: {ucr/mon:.2f},   USP speedup: {usp/mon:.2f}\n")


    # Compare timing: when is UCR MON slower than UCR and UCR USP ?
    def slower_than(name1, name2, doPrint=False):
        # Count total number of analysed result
        count_total = 0
        # count name1>name2
        count = 0
        # diff name1-name
        diff_total = 0.0
        diff_min = float("inf")
        diff_max = 0.0
        # runtimes
        n1total = 0.0
        n2total = 0.0
        # List of datasets where name1>name2
        slowerlist = []

        # Iterate over datasets
        for d in datasets:
            for q in queries:
                for l in lengths:
                    for w in ratios:
                        count_total += 1
                        n1 = data[name1][d][q][l][w][0]
                        n2 = data[name2][d][q][l][w][0]
                        n1total += n1
                        n2total += n2
                        if n1 > n2:
                            count+=1
                            diff = n1 - n2
                            diff_total+=diff
                            diff_min = min(diff_min, diff)
                            diff_max = max(diff_max, diff)
                            slowerlist.append((d, q, l, w, n1, n2))
        # Print slower cases
        if doPrint and slowerlist:
            print(name1 + " slower than " + name2 + ":")
            for (d, q, l, w, n1, n2) in slowerlist:
                diff = n1-n2
                ratio = (100*diff)/n2
                head = (str(d) + " " + str(q) + " " +str(l) + " " +str(w))
                print(f"  {head:<25}  {n1:8.3f} > {n2:<8.3f}   diff: {diff:<8.3f}  slower + %: {ratio:.3f}")

        # Return info
        return (count_total, count, diff_total, diff_min, diff_max, name1, name2, n1total, n2total)


    # Printing infor return by the above
    def print_info(info):
        (count_total, count, diff, diff_min, diff_max, name1, name2, n1total, n2total) = info
        print(f"Number of datasets where {name1} is slower than {name2}: {count} / {count_total} =  {float(count/count_total):5f}")
        print(f"  Average slower by: {diff/count:.5f}s")
        print(f"  Min slower by:     {diff_min:8.5f}s")
        print(f"  Max slower by:     {diff_max:8.5f}s")
        print("Total runtimes:")
        print(f"  {name1:<10}: {n1total:.2f}s = {datetime.timedelta(seconds=n1total)}")
        print(f"  {name2:<10}: {n2total:.2f}s = {datetime.timedelta(seconds=n2total)}")
        txt = name1+"/"+name2
        print(f"  {txt:<10}: {n1total/n2total:.3f}")
        txt = name2+"/"+name1
        print(f"  {txt:<10}: {n2total/n1total:.3f}")


    print("--------------------------------------------")
    print("Comparison: who is slower than the other?")

    doPrint=False

    # MON > UCR
    print()
    print_info(slower_than("UCR_MON", "UCR", doPrint))

    # MON > USP
    print()
    print_info(slower_than("UCR_MON", "UCR_USP", doPrint))

    print()
    print()

    # MON_nolb > UCR
    print()
    print_info(slower_than("UCR_MON_nolb", "UCR", doPrint))

    # MON_nolb > USP
    print()
    print_info(slower_than("UCR_MON_nolb", "UCR_USP", doPrint))

    print()
    print()

    # USP > UCR
    print()
    print_info(slower_than("UCR_USP", "UCR", doPrint))
