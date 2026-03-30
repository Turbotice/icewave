import matplotlib.pyplot as plt

# fonction pour mesurer les coordonnées de n points sur une figure

def get_n_points(x, y, n_points=1):
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.grid(True)

    coords = []

    def onclick(event):
        # 🚫 Ignore si zoom ou pan actif
        if fig.canvas.toolbar.mode != '':
            return

        # 🚫 Ignore clic hors graphe
        if event.xdata is None or event.ydata is None:
            return

        # 🚫 Seulement clic gauche
        if event.button != 1:
            return

        coords.append((event.xdata, event.ydata))
        print(f"[{len(coords)}/{n_points}] x = {event.xdata:.4f}, y = {event.ydata:.4f}")

        ax.plot(event.xdata, event.ydata, 'ro')
        fig.canvas.draw()

        if len(coords) >= n_points:
            print("✅ Terminé")
            fig.canvas.mpl_disconnect(cid)
            plt.close(fig)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    print(f"👉 Clique sur {n_points} point(s)")
    print("👉 Utilise zoom/pan SANS déclencher de clic parasite")

    plt.show()
    return coords