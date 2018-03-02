import os
from math import *
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
from math import log
import argparse


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
    return c * r * 1


def get_kl_dis(dict1, dict2):
    dis = 0
    for k in dict1:
        dis += dict1[k] * log(dict1[k] * 1.0 / dict2[k])
    return -dis


def lrcf_weight(loc_trace, isp_traces):

    loc_points = dict(zip([k.split(",")[0] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0],
                          [k.split(",")[1] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0]))
    loc_region_dict = dict()
    for idx in loc_points:
        min_dis = 1.001
        min_region = ''
        loc_loc = loc_points[idx]
        for region in region_c_dict:
            loc = region_c_dict[region]
            dis = get_distance(float(loc.split("_")[0]), float(loc.split("_")[1]), float(loc_loc.split("_")[1]),
                               float(loc_loc.split("_")[0]))
            if dis < min_dis:
                min_dis = dis
                min_region = region
        if min_region in region_c_dict:
            if idx in loc_region_dict:
                loc_region_dict[idx].append(region_c_dict[min_region])
            else:
                loc_region_dict[idx] = [region_c_dict[min_region]]
    weight_list = []
    isp_ids = []
    loc_region_vec = dict()
    for region in loc_region_dict:
        region_mapped = loc_region_dict[region]
        for region_mapped_one in region_mapped:
            if region_mapped_one in loc_region_vec:
                loc_region_vec[region_mapped_one] += 1
            else:
                loc_region_vec[region_mapped_one] = 1

    for region in loc_region_vec:
        if region in region_mapped_count:
            # print 'find'
            loc_region_vec[region] = loc_region_vec[region] * 1.0 / (1.0 + log(region_mapped_count[region]))
        else:
            pass

    loc_regiones = loc_region_vec.keys()

    for isp_trace in isp_traces:
   
        isp_id = isp_trace.split(";")[0]
        isp_ids.append(isp_id)
        try:
            isp_points = dict(zip([k.strip().split(",")[0] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2],
                                  [k.strip().split(",")[1] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2]))
            isp_region_vec = dict()
            for region in isp_points:
                if isp_points[region] not in centre_map_dict:
                    continue
                region_mapped = centre_map_dict[isp_points[region]]
                if region_mapped in isp_region_vec:
                    isp_region_vec[region_mapped] += 1
                else:
                    isp_region_vec[region_mapped] = 1
           
            # isp_region_vec_valuesum = sum(isp_region_vec.values()) * 1.0
            
            # for k in isp_region_vec:
            #     isp_region_vec[k] /= isp_region_vec_valuesum
            for region in isp_region_vec:
                if region in region_mapped_count:
                    # print 'in region_count_dict'
                    isp_region_vec[region] = isp_region_vec[region] * 1.0 / ( 1.0 + log(region_mapped_count[region]))
                else:
                    print 'not in region_count_dict'
            
            
            common_region = filter(lambda x:x in isp_region_vec, loc_regiones)
            if len(common_region) > 0:
                for region in loc_region_vec:
                    if region not in isp_region_vec:
                        isp_region_vec[region] = 0.0
                for region in isp_region_vec:
                    if region not in loc_region_vec:
                        loc_region_vec[region] = 0.0

                loc_region_vec_list = [loc_region_vec[k] for k in loc_region_vec]
                isp_region_vec_list = [isp_region_vec[k] for k in loc_region_vec]

                if len(loc_region_vec_list) == 0 or len(isp_region_vec_list) == 0 :
                    cos_dis = -1.0
                else:
                    cos_dis = cosine_similarity(np.array(loc_region_vec_list).reshape(1,-1), np.array(isp_region_vec_list).reshape(1,-1)).sum()

            else:
                cos_dis = -1.0
            weight_list.append((isp_id, cos_dis))
        except IndexError:
            print 'error', isp_trace

    return weight_list


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path,save_weights=False):
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0]
    isp_traces = [all_traces[k] for k in range(len(all_traces)) if k >= 1]
    ispUID_score = lrcf_weight(loc_trace, isp_traces)
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
    parser.add_argument('--output_path', nargs='?', default='./output/lrcf/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters',
                        help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=10, help='K in hit-precision', type=int)
    parser.add_argument('--region_popularity', nargs='?', default='./RegionPopularity', help='Popularity of difference regions')

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

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)