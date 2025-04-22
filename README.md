# astro_comp_group
**Proyecto 03**

**Integrantes:**
- David Santiago Cuchigay 
- Daniel Rojas
- David Sáenz 
- Óscar Calvo
## Respuestas preguntas
**3.¿El comportamiento mejora o empeora al modificar el paso de la integración?**

— El comportamiento mejora siempre y cuando el tamaño del paso se encuentre a una distancia confiable del cero de la máquina, pues al esta reducirse considerablemente, el algoritmo incrementará su error debido a redondeo, ruido numérico y acumulación de los mismos.

**5.¿Mejora el tiempo de computo utilizando este método? ¿Se conserva la energía total del sistema?**

Usando los datos en SOstars.dat, se simulan cien años con 100000 pasos, utilizando los siguientes métodos: Runge-Kutta Orden 4 (RK4), Velocity Verlet y Bulirsch-Stoer. Los tiempos de integración, y los cambios de energía obtenidos son:

## Resultados
### Tiempo de integración
| Método              | Tiempo de integración |
|---------------------|----------------------|
| RK4                 | 143.045 s            |
| Velocity Verlet     | 89.677 s             |
| Bulirsch-Stoer      | 16.471 s             |

### Conservación de la energía
| Método              | Cambio de energía $\Delta E$   |
|---------------------|----------------------|
|RK4                  | 0.67535              |
|Velocity Verlet      | -0.0388              |
|Bulirsch-Stoer       | -0.0001148           |


El método de Bulirsch-Stoer es más rápido y conserva mejor la energía total del sistema, ya que tiene un error de energía dos y tres ordenes de magnitud menor que los métodos de Runge Kutta 4 y Velocity Verlet, respectivamente.
