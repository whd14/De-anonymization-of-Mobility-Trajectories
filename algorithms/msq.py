import os
from math import *
import argparse


def get_distance(lon1, lat1, lon2, lat2):
    # return distance (kilometer) of two GPS locations
    lon1 = float(lon1)
    lat1 = float(lat1)
    lon2 = float(lon2)
    lat2 = float(lat2)
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6378.137
    return c * r 


def msq_weight(loc_trace, isp_traces):
    Dmax=20
    try:
        loc_points = dict(zip([k.strip().split(",")[0] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0],
                              [k.strip().split(",")[1] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0]))
    except IndexError:
        print 'loc_trace', loc_trace
        print 'isp_trace', isp_trace


    weight_list = []
    for isp_trace in isp_traces:
        isp_id = isp_trace.split(";")[0]
        try:
            isp_points = dict(zip([k.strip().split(",")[0] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2],
                                  [k.strip().split(",")[1] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2]))
            weight=0
            for idx in loc_points.keys():
                locLOC=loc_points[idx]
                if isp_points.has_key(idx):
                    rg=isp_points[idx]
                    locRG = region_c_dict[rg]
                    dis = get_distance(float(locLOC.split("_")[1]), float(locLOC.split("_")[0]), float(locRG.split("_")[0]),
                                       float(locRG.split("_")[1]))
                    if dis>=Dmax:
                        dis=Dmax
                    weight=weight-dis**2
                else:
                    weight=weight-Dmax**2
            weight_list.append((isp_id,weight))
        except IndexError:
            print 'data error'
    return weight_list


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path, save_weights=False):
    # deanonymization task of single user
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0]
    isp_traces = [all_traces[k] for k in range(len(all_traces)) if k >= 1]
    ispUID_score = msq_weight(loc_trace, isp_traces)
    ispUID_score_sorted = sorted(ispUID_score, key=lambda x: x[1], reverse=True)
    ispUID_sorted = [x[0] for x in ispUID_score_sorted]
    locUID = loc_trace.split(";")[0]
    if locUID in ispUID_sorted:
        rank_true = ispUID_sorted.index(locUID)
    else:
        rank_true = len(ispUID_sorted)
    acc = hit_precision(rank_true, K_hp)
    if save_weights:
        fw = open(output_path + locfile_name, 'w')
        for i in range(len(ispUID_score_sorted)):
            fw.write(str(ispUID_score_sorted[i][0]) + '\t' + str(ispUID_score_sorted[i][1]) + '\n')
        fw.close()
    return locfile_name, acc


def main(input_path='./input/', output_path='./output/', save_weights=False):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    file_names = os.listdir(input_path)
    acc_list = []
    for file_name in file_names:
        acc_list.append(single_task(file_name, input_path, output_path, save_weights))
    with open(output_path + '/result', 'w') as fw:
        for acc in acc_list:
            fw.write(str(acc[0]) + '\t' + str(acc[1]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', nargs='?', default='./input/', help='Path of input data')
    parser.add_argument('--output_path', nargs='?', default='./output/msq/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters', help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=10, help='K in hit-precision', type=int)

    args = parser.parse_args()
    K_hp = args.K
    region_centre = args.region_centre
    region_c_list = open(region_centre, 'r').readlines()
    region_c_dict = dict(
        zip([k.split("_")[0] for k in region_c_list], [k.split("_")[1] + "_" + k.split("_")[2] for k in region_c_list]))

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)
