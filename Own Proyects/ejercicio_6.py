# Ejercicio 6: LISTAS
# Autor: Alejandro Antonio Espinosa Gil

# Escribir un programa que almacene las asignaturas de un curso (por ejemplo Matemáticas,
# Física, Química, Historia y Lengua) en una lista, pregunte al usuario la nota que ha
# sacado en cada asignatura y elimine de la lista las asignaturas aprobadas.
# Al final el programa debe mostrar por pantalla las asignaturas que el usuario tiene
# que repetir.


# Lista de asignaturas
matList = list(("Matemáticas","Física","Química","Historia","Lengua"))

# Lista donde se almacenarán las asignaturas aprobadas y que serán eliminadas
# del la lista inicial de asignaturas
matAp = []

# Ciclo for que abarca el tamaño de la lista de asignaturas
for i in range(0,len(matList),1):
    # El usuario introduce la nota de la materia correspondiente
    calif = int(input("Ingresar la nota de "+matList[i]+": "))

    # Ciclo while que obliga al usuario a introducir un número válido de calificación
    # entre 0 y 100
    while calif > 100 or calif < 0:
        print("Introduzca un valor de la calificación entre ")
        calif = int(input("Ingresar la nota de "+i+": "))
    
    # La calificación aprobatoria es 60, por lo tanto, si la calificación es mayor
    # a este valor, se considera aprobada la asignatura
    if calif > 59:
        # Se almacenan las asignaturas aprobadas en una lista
        matAp.append(matList[i])

# Filtro que sirve para eliminar las asignaturas aprobadas de la lista
# original de asignaturas
matRep = [i for i in matList if i not in matAp]


if len(matRep) == 0: # Si no hay ninguna asignatura en la lista notifica al estudiante
    print("No tienes asignaturas por repetir.")

else: # Se imprimen las asignaturas que se deben repetir
    print("Las asignaturas por repetir son: ")
    for i in matRep:
        print("- "+i)


# Esta es una prueba para ver si se actualiza en github