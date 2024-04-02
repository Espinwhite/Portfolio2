# Ejercicio 12 Diccionario
# Autor: Alejandro Antonio Espinosa Gil

# Escribir un programa que almacene el diccionario con los créditos de las asignaturas
# de un curso {'Matemáticas': 6, 'Física': 4, 'Química': 5} y después muestre por pantalla
# los créditos de cada asignatura en el formato <asignatura> tiene <créditos> créditos,
# donde <asignatura> es cada una de las asignaturas del curso, y <créditos> son sus créditos.
# Al final debe mostrar también el número total de créditos del curso.

# Diccionario de asignaturas con créditos
asignaturas = {'Matematicas': 6, "Física": 4, "Química": 5}
totcred = 0             # Contador de créditos totales
for i in asignaturas:
    mat = asignaturas[i]    # Se obtiene el nombre de la materia
    cred = asignaturas.get(i)   # Se obtiene el crédito correspondiente a cada materia
    print(i +" tiene "+str(cred)+" créditos")
    totcred = totcred + cred    # Se suman los créditos de cada materia en cada iteración

print("Los créditos totales del curso son: ", totcred)