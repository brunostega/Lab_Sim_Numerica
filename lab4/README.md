Dinamica molecolare, Algoritmo di Verlet


Per compilare ed eseguire il programma per tutte e tre le fasi con equilibrazione, lanciare il comando
>./experiment_with_equilibration.sh

Se il sistema è già stato equilibrato (vedere file "output.temp" nelle tre cartelle), è possibile evitare la termalizzazione usando il comando

>./experiment.sh

In questo caso il programma partirà dalla configurazione opportunamente equilibrata per le diverse fasi.

Nel caso in cui si volesse eseguire il programa per una singola fase (gas, liquido, solido), è possibile lanciare i comandi

>make Gas

>make Liquid

>make Solid

nel caso in cui si volesse equilibrare il sistema, e

>make Gas_termalized

>make Liquid_termalized

>make Solid_termalized


