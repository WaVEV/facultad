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



Toma del pulso de mejoras:
	Compilando con el flag de optimizacion O3 y en "mini", con los siguientes parametros:

		EPS= 3.0e-14
		R_MIN= 0.0 
		R_MAX= 50.0
		L_INTERVALS= 510
		KORD= 5 
		RADIO_1= 5.0 
		RADIO_2= 10.0 
		ME= 1.0 
		INT_G= 500 
		NEV= 15 
		L_MAX= 0 
		LAMBDA_IN= 0.0 
		LAMBDA_FIN= 20.0 
		NUMEROS_PUNTO_LAMBDA= 200 
		BASE_KORD= 0

	Performance:
		Sólo se toma el tiempo que tarda el proceso calculo_matrices
		Código base:                            2.727035049000 segundos
		Con cacheo de la funcion bsplvb: 	2.177257117000 segundos
                Reemplazando los arrelgos t y k:        1.811051266995 segundos
                Recalculo de un indice:                 1.772214113967 segundos
        

Apéndice:
        Lo primero que se hará para intentar paralelizar es usar algo asi como una "ventana corredisa" (en verdad esto fue lo primero que se me ocurrió) para operar en el circulo KNOTS para que no haya condicion de carrera cuando varios cores intenten escribir en el mismo lugar.


PD: le pasaré un corrector, estoy haciendo esto mas para mi que para presentar, y tambien se que
que le sirve al lector.
