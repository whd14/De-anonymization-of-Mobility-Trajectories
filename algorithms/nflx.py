# -*- coding:utf-8 -*-
import os
from math import *
import argparse

#=============================hyper-parameter===========================
rou0 = 10
n0 = 2
#=======================================================================


def get_distance(lon1, lat1, lon2, lat2):
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


def nflx_weight(loc_trace, isp_traces):
    weight = []
    loc_trace_dict = dict()
    id_points = loc_trace.split(";")
    points = id_points[1].split("|")
    for point in points:
        if len(point.strip()) > 0:
            time_loc = point.split(",")
            latlng = time_loc[1]
            min_dis = 1.001
            min_region = ''
            for base_one in region_c_dict:
                if base_one not in centre_map_dict:
                    continue
                dis = get_distance(latlng.split("_")[1], latlng.split("_")[0],
                                   region_c_dict[base_one].split("_")[0], region_c_dict[base_one].split("_")[1])
                if dis < min_dis:
                    min_dis = dis
                    min_region = base_one
            if min_region in centre_map_dict:
                region = centre_map_dict[min_region]
                if region in loc_trace_dict:
                    loc_trace_dict[region].append(time_loc[0])
                else:
                    loc_trace_dict[region] = [time_loc[0]]
    # if len(loc_trace_dict) == 0:
    #     return locfile_name, 0
    weight_list = []
    for isp_trace in isp_traces:
        isp_trace_dict = dict()
        id_points = isp_trace.split(";")
        id = id_points[0]
        points = id_points[1].split("|")
        for point in points:
            if len(point.strip()) > 0:
                time_region_loc = point.split(",")
                if len(time_region_loc) < 2:
                    continue
                if time_region_loc[1] not in centre_map_dict:
                    continue
                region = centre_map_dict[time_region_loc[1]]
                if region in isp_trace_dict:
                    isp_trace_dict[region].append(time_region_loc[0])
                else:
                    isp_trace_dict[region] = [time_region_loc[0]]
        trace1 = loc_trace_dict
        trace2 = isp_trace_dict
        weight = 0
        for region in trace1:
            timelist_1 = trace1[region]
            if region in trace2:
                average_dis = 0
                timelist_2 = trace2[region]
                for time in timelist_1:
                    timedis = [abs(float(k) - float(time)) for k in timelist_2]
                    try:
                        mindis = min(timedis)
                        average_dis += mindis
                    except ValueError:
                        pass
                average_dis /= len(timelist_1)
                weight_one = exp(len(timelist_1) / n0) + exp(- average_dis / len(timelist_1) / rou0)
                weight += weight_one
        weight_list.append((id, weight))
    return weight_list


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
    ispUID_score = nflx_weight(loc_trace, isp_traces)
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
    parser.add_argument('--output_path', nargs='?', default='./output/nflx/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters',
                        help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=100000, help='K in hit-precision', type=int)

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
