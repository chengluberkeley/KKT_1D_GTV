#echo "=== Start to run l2-l1-nw ==="
./kkt_main -profile l2-l1-nw
#echo "=== Complete l2-l1-nw ==="

#echo "=== Start to run l2-l1-w ==="
./kkt_main -profile l2-l1-w
#echo "=== Complete l2-l1-w ==="

#echo "=== Start to run l2-l2 ==="
./kkt_main -profile l2-l2
#echo "=== Complete l2-l2 ==="

#echo "=== Start to run lp-lq ==="
./kkt_main -profile lp-lq > lp_lq_runtimes.txt
#echo "=== Complete lp-lq ==="

#echo "=== Start to run pwl2-l1 ==="
./kkt_main -profile pwl2-l1 -num_scales 6
#echo "=== Complete pwl2-l1 ==="

#echo "=== Start to run pwl1-l1 ==="
./kkt_main -profile pwl1-l1 -nloglogn_stop_scale 4 -fix_n 10000 -num_scales 6
#echo "=== Complete pwl1-l1 ==="

#echo "=== Start to run l1-l1 ==="
./kkt_main -profile l1-l1 -nloglogn_stop_scale 4 -fix_n 10000
#echo "=== Complete l1-l1 ==="

echo "=== Start to run Huber ==="
./kkt_main -profile Huber
echo "=== Complete Huber ==="

