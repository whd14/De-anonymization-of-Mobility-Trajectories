# -*- coding:utf-8 -*-
import os
from math import *
import argparse

# =============================hyper-parameter===========================
p1 = 0.3
p2 = 0.01
totalT = 192.0
totalN = 1798550.0
# =======================================================================


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


def weight_00(p):
    return log((1 - p + p * (1 - p1) * (1 - p2)) / ((1 - p + p * (1 - p1)) * (1 - p + p * (1 - p2))))


def weight_01(p):
    return log((1 - p1) / (1 - p * p1)) - weight_00(p)


def weight_10(p):
    return log((1 - p2) / (1 - p * p2)) - weight_00(p)


def weight_11(p):
    return log(1 / p) - weight_00(p)


def add(x, y):
    return x + y


def pois_weight(loc_trace, isp_traces):
    loc_rec = loc_trace.split(";")[1]
    loc_idx_latlng_list = loc_rec.split("|")
    loc_idx_point = dict(zip([k.split(",")[0] for k in loc_idx_latlng_list if len(k) > 0],
                             [k.split(",")[1] for k in loc_idx_latlng_list if len(k) > 0]))
    loc_idx_region = dict()
    for loc_idx in loc_idx_point:
        latlng = loc_idx_point[loc_idx]
        lat = latlng.split('_')[0]
        lng = latlng.split('_')[1]
        min_dis = 1.001
        min_region = ''
        for region in region_c_dict:
            region_lat = region_c_dict[region].split('_')[1]
            region_lng = region_c_dict[region].split('_')[0]
            dis = get_distance(region_lng, region_lat, lng, lat)
            if dis <= min_dis:
                min_dis = dis
                min_region = region
        if len(min_region) == 0:
            continue
        else:
            loc_idx_region[loc_idx] = centre_map_dict[min_region]
    isp_recs = [k.split(";")[1] for k in isp_traces]
    weight_list = []
    for isp_rec in isp_traces:
        isp_id = isp_rec.split(';')[0]
        isp_rec = isp_rec.split(';')[1]
        weight = 0
        coincidence = dict()
        isp_idx_region_latlng_list = isp_rec.split("|")
        isp_idx_point = dict(
            zip([k.split(",")[0] for k in isp_idx_region_latlng_list if len(k) > 0 and len(k.split(",")) == 2],
                [k.split(",")[1] for k in isp_idx_region_latlng_list if len(k) > 0 and len(k.split(",")) == 2]))
        isp_idx_region = dict()
        for idx in isp_idx_point:
            region = isp_idx_point[idx]
            if region in centre_map_dict:
                isp_idx_region[idx] = centre_map_dict[region]
        for idx in loc_idx_region:
            loc_region = loc_idx_region[idx]
            if idx in isp_idx_region:
                if loc_region == isp_idx_region[idx]:
                    coincidence[idx] = loc_region
                    weight += weight_11(region_p_dict[coincidence[idx]])

        isp_none = [weight_10(region_p_dict[isp_idx_point[k]]) for k in isp_idx_point if
                    k not in coincidence
                    and region_p_dict.has_key(isp_idx_point[k])]
        weight += sum(isp_none)
        loc_none = [weight_01(region_p_dict[loc_idx_point[k]]) for k in loc_idx_point if k not in coincidence
                    and region_p_dict.has_key(loc_idx_point[k])]
        weight += sum(loc_none)
        weight_list.append((isp_id, weight))
    return weight_list


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path, save_weights=False):
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0].strip()
    isp_traces = [all_traces[k].strip() for k in range(len(all_traces)) if k >= 1]
    ispUID_score = pois_weight(loc_trace, isp_traces)
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
    parser.add_argument('--output_path', nargs='?', default='./output/hist/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters',
                        help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=10, help='K in hit-precision', type=int)
    parser.add_argument('--region_popularity', nargs='?', default='./RegionPopularity',
                        help='Popularity of difference regions')

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

    region_pop_list = open(args.region_popularity, 'r').readlines()
    region_pop_dict = dict(
        zip([k.strip().split()[0] for k in region_pop_list], [k.strip().split()[1] for k in region_pop_list]))
    region_mapped_count = dict()
    for region in centre_map_dict:
        if region not in region_pop_dict:
            continue
        region_mapped = centre_map_dict[region]
        if region_mapped in region_mapped_count:
            region_mapped_count[region_mapped] += int(region_pop_dict[region])
        else:
            region_mapped_count[region_mapped] = int(region_pop_dict[region])
    region_p_dict = region_mapped_count
    for k in region_p_dict:
        region_p_dict[k] = region_p_dict[k] * 1.0 / (totalN * totalT)
    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)
