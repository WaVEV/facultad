Aca iré escribiendo lo que hago día a día

23/06/2015

    Tareas hechas:
        1) Se cambiaron muchos nombres de las variables globales
        2) Se definió como estatica el arreglo Sp
        3) Se hicieron variables globales con los arreglos que se definian en el main

    Problemas:
        1) Los arreglos del punto 3 de tarea se pasaron a global porque los 
            declaraba en el main y se le daba memoria en tiempo de ejecución, 
            y es mas, daba sigsfaul. tenes alguna idea de que puede ser

    Tareas para la proxima:
        1) Hacer que las variables globales sean globales y no se esten pasando por argunemnto
        2) Hacer el dibujo de como opera la matriz
        3) yapa: perf

28/06/2015
    Tareas hechas:
        1) Se calculo de manera O(1) (con operaciones simpples) los arreglos t y k
        2) Se realizo el pasaje (o no pasaje) de parametros globales
        3) Se empezo a hacer el dibujo de la operacion matricial

    Problemas:
        --

    Tareas para la próxima:
        1) Terminar el dibujo de la addicion en la matriz
        2) Implemetar la primera idea de paralelización
        3) Usar los pragmas OPM






Apéndice:
        Lo primero que se hará para intentar paralelizar es usar algo asi como una "ventana corredisa" (en verdad esto fue lo primero que se me ocurrió) para operar en el circulo KNOTS para que no haya condicion de carrera cuando varios cores intenten escribir en el mismo lugar.


PD: le pasaré un corrector, estoy haciendo esto mas para mi que para presentar, y tambien se que
que le sirve al lector.
