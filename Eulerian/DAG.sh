euler_library="/home/qlalf1457/InfoGenomeR/Eulerian/DAG_entropy"

parameters=`cat edge_information.txt   | awk 'BEGIN{total_nodes=0;total_edges=0; hidden_node=0; hidden_edge=0}{if($1>total_nodes){total_nodes=$1};if($2>total_nodes){total_nodes=$2};total_edges=$6;hidden_edge=$7}END{print (total_nodes+2)"\t"(total_edges+1)"\t"(total_nodes+1)"\t"(hidden_edge+1)}'`
total_nodes=`echo $parameters |awk '{print $1}'`;
total_edges=`echo $parameters| awk '{print $2}'`;
hidden_node=`echo $parameters| awk '{print $3}'`;
hidden_edge=`echo $parameters| awk '{print $4}'`;

if [[ $1 == "T" ]];then
	level_index=`$euler_library/main true 0 0 $total_nodes $total_edges $hidden_node $hidden_edge | grep level_index`;
	index1=`echo $level_index | awk '{print $2}'`;
	index2=`echo $level_index | awk '{print $3}'`;
	
	i=$index1;
	while [[ $i -lt $index2 ]];do
		$euler_library\/main false $i $index2 $total_nodes $total_edges $hidden_node $hidden_edge | grep new;
		i=$(($i+1));
	done
else

	$euler_library/main false 0 0 $total_nodes $total_edges $hidden_node $hidden_edge 
fi

