# astro_comp_group
**Proyecto 03**

**Integrantes:**
- David Santiago Cuchigay 
- Daniel Rojas
- David Sáenz 
- Óscar Calvo
## Respuestas preguntas
**3.¿El comportamiento mejora o empeora al modificar el paso de la integración?**

—Esto depende, porque aunque al principio el error se reduce cuando se incrementa el número de pasos, llega un momento en que esto se invierte. Lo anterior ocurre por el "cero de la máquina". A partir de esto, si se siguen aumentando los pasos, el error empieza a crecer nuevamente debido al ruido numérico o a los errores de redondeo que se acumulan por las limitaciones de la máquina.

**5.¿Mejora el tiempo de computo utilizando este método? ¿Se conserva la energía total del sistema?**

Usando los datos en SOstars.dat, se simularon cien años con 100000 pasos, utilizando los siguientes métodos: Runge-Kutta Orden 4 (RK4), Velocity Verlet y Bulirsch-Stoer. Los tiempos de integración, y los cambios de energía fueron:

## Resultados
### Tiempo de integración
| Método              | Tiempo de integración |
|---------------------|----------------------|
| RK4                 | 143.04559588 s       |
| Velocity Verlet     | 89.67779040336609 s  |
| Bulirsch-Stoer      | 16.4712 s            |

### Conservación de la energía
| Método              | Cambio de energía $\delta E$   |
|---------------------|----------------------|
|RK4                  | 0.67535              |
|Velocity Verlet      | -0.0388              |
|Bulirsch-Stoer       | -0.0001148           |


El método de Bulirsch-Stoer es más rápido y conserva mejor la energía total del sistema, ya que tiene un error de energía mucho menor que los otros dos métodos.
