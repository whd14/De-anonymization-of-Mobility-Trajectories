# -*- coding:utf-8 -*-
import os
from math import *
import argparse


def get_distance(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6378.137
    return c * r


def wyci(loc_trace, isp_traces):
    alpha = 0.1
    probability = []
    loc_points = [k.split(",")[1] for k in loc_trace.strip().split(";")[1].split("|") if len(k) > 0]
    loc_region = []
    for loc_point in loc_points:
        min_dis = 1.001
        min_region = ''
        for region in region_c_dict:
            loc = region_c_dict[region]
            dis = get_distance(float(loc.split("_")[0]), float(loc.split("_")[1]), float(loc_point.split("_")[1]), float(loc_point.split("_")[0]))
            if dis < min_dis:
                min_dis = dis
                min_region = region
        if min_region in centre_map_dict:
            loc_region.append(centre_map_dict[min_region])

    for isp_trace in isp_traces:
        isp_id = isp_trace.split(";")[0]
        isp_points = isp_trace.strip().split(";")[1].split("|")
        points_num = len(isp_points)
        isp_freq = dict()
        for isp_point in isp_points:
            if len(isp_point) > 0:
                idx_region_loc = isp_point.split(",")
                if len(idx_region_loc) < 2:
                    continue
                if idx_region_loc[1] not in centre_map_dict:
                    continue
                if centre_map_dict[idx_region_loc[1]] in isp_freq:
                    isp_freq[centre_map_dict[idx_region_loc[1]]] += 1.0
                else:
                    isp_freq[centre_map_dict[idx_region_loc[1]]] = 1.0
        min_point = '-1'
        p = 1
        for loc_region_one in loc_region:
            if loc_region_one in isp_freq:
                p *= (isp_freq[loc_region_one] + alpha) / (points_num + alpha * 3000)
            else:
                p *= (0 + alpha) / (points_num + alpha * 3000)
        probability.append((isp_id, p))
    return probability


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path, save_weights=False):
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0]
    isp_traces = [all_traces[k] for k in range(len(all_traces)) if k >= 1]
    ispUID_score = wyci(loc_trace, isp_traces)
    ispUID_score_sorted = sorted(ispUID_score, key=lambda x: x[1], reverse=True)
    ispUID_sorted = [x[0] for x in ispUID_score_sorted]
    locUID = loc_trace.split(";")[0]
    if locUID in ispUID_sorted:
        rank_true = ispUID_sorted.index(locUID)
    else:
        rank_true = len(ispUID_sorted)
    print rank_true
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
    parser.add_argument('--output_path', nargs='?', default='./output/wyci/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters',
                        help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=100, help='K in hit-precision', type=int)

    args = parser.parse_args()
    K_hp = args.K
    region_centre = args.region_centre
    region_c_list = open(region_centre, 'r').readlines()
    region_c_dict = dict(
        zip([k.split("_")[0] for k in region_c_list], [k.split("_")[1] + "_" + k.split("_")[2] for k in region_c_list]))

    centre_map_dict = {}
    for centre in region_c_dict:
        loc = region_c_dict[centre].split('_')
        new_loc = str(floor(float(loc[0]) / args.scale) * args.scale) + '-' + str(
            floor(float(loc[1]) / args.scale) * args.scale)
        centre_map_dict[centre] = new_loc

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)
