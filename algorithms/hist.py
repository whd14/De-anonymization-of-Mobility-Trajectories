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


def get_kl_dis(dict1, dict2):
    # return KL-distance of two vectors
    dis = 0
    for k in dict1:
        dis += dict1[k] * log(dict1[k] * 1.0 / dict2[k])
    return -dis


def hist_weight(loc_trace, isp_traces):
    # return HIST weight of a LOC-trace and several ISP-traces
    loc_points = dict(zip([k.split(",")[0] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0],
                          [k.split(",")[1] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0]))
    loc_centre_dict = dict()
    for idx in loc_points:
        min_dis = 1.001
        min_centre = ''
        loc_loc = loc_points[idx]
        for centre in region_c_dict:
            loc = region_c_dict[centre]
            dis = get_distance(float(loc.split("_")[0]), float(loc.split("_")[1]), float(loc_loc.split("_")[1]),
                               float(loc_loc.split("_")[0]))
            if dis < min_dis:
                min_dis = dis
                min_centre = centre
        if min_centre in centre_map_dict:
            if idx in loc_centre_dict:
                loc_centre_dict[idx].append(centre_map_dict[min_centre])
            else:
                loc_centre_dict[idx] = [centre_map_dict[min_centre]]
    weight_list = []
    loc_centre_vec = dict()
    for centre in loc_centre_dict:
        centre_mapped = loc_centre_dict[centre]
        for centre_mapped_one in centre_mapped:
            if centre_mapped_one in loc_centre_vec:
                loc_centre_vec[centre_mapped_one] += 1
            else:
                loc_centre_vec[centre_mapped_one] = 1
    loc_centre_vec_valuesum = sum(loc_centre_vec.values()) * 1.0
    for k in loc_centre_vec:
        loc_centre_vec[k] /= loc_centre_vec_valuesum
    for isp_trace in isp_traces:
        isp_id = isp_trace.split(";")[0]
        try:
            isp_points = dict(zip([k.strip().split(",")[0] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2],
                                  [k.strip().split(",")[1] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2]))
            isp_centre_vec = dict()
            for centre in isp_points:
                if isp_points[centre] not in centre_map_dict:
                    continue
                centre_mapped = centre_map_dict[isp_points[centre]]
                if centre_mapped in isp_centre_vec:
                    isp_centre_vec[centre_mapped] += 1
                else:
                    isp_centre_vec[centre_mapped] = 1

            isp_centre_vec_valuesum = sum(isp_centre_vec.values()) * 1.0
            for k in isp_centre_vec:
                isp_centre_vec[k] /= isp_centre_vec_valuesum
            merge_dict = dict(isp_centre_vec)
            for loc_i in loc_centre_vec:
                if loc_i in merge_dict:
                    merge_dict[loc_i] += loc_centre_vec[loc_i]
                else:
                    merge_dict[loc_i] = loc_centre_vec[loc_i]
            for merge_i in merge_dict:
                merge_dict[merge_i] *= 0.5
            kl_dis = get_kl_dis(isp_centre_vec, merge_dict) + get_kl_dis(loc_centre_vec, merge_dict)
            weight_list.append((isp_id,kl_dis))
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
    ispUID_score = hist_weight(loc_trace, isp_traces)
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
    parser.add_argument('--output_path', nargs='?', default='./output/hist/', help='Path of output data')
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

    centre_map_dict = {}
    for centre in region_c_dict:
        loc = region_c_dict[centre].split('_')
        new_loc = str(floor(float(loc[0]) / args.scale) * args.scale) + '-' + str(floor(float(loc[1]) / args.scale) * args.scale)
        centre_map_dict[centre] = new_loc

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)
