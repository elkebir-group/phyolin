 
#!/bin/bash
columns=`printf "%s\t%s" "file seconds"`

echo $columns > ../runtimes.txt


 for f in ./*.csv; do
          echo $f
          fname=$(echo "$f" | grep -oP 'AML-[0-9]*_rep[0-9]*.csv')
         
            ofile="../output/${fname}"
        
            echo $ofile
       
   
  

        #run each file and save the run times
        fstart=$(date +%s.%N)

            ../../../build/phyolin $f $ofile
         
        phyduration=$(echo "$(date +%s.%N) - $fstart" | bc)
        exec_time=`printf "%s\t%.2f" $f $phyduration`
        echo $exec_time >> ../runtimes.txt
    done 