def cut_link_in_mesh(filename_in, filename_out, v1, v2):
    """
    Supprime tous les éléments (Hexahedra, Tetrahedra, ...) du fichier .mesh
    qui contiennent à la fois les sommets v1 et v2.

    Parameters
    ----------
    filename_in : str
        Fichier .mesh d'entrée.
    filename_out : str
        Nouveau fichier .mesh de sortie.
    v1, v2 : int
        Indices des sommets à délier.
    """
    with open(filename_in, "r") as f:
        lines = f.readlines()

    out_lines = []
    in_block = None
    n_removed = 0

    for line in lines:
        # Détection début/fin des blocs
        if line.strip().startswith("Hexahedra") or line.strip().startswith("Tetrahedra"):
            in_block = line.split()[0]  # "Hexahedra" ou "Tetrahedra"
            out_lines.append(line)
            continue
        elif line.strip() == "End":
            in_block = None
            out_lines.append(line)
            continue

        # Si on est dans un bloc de connectivités
        if in_block:
            parts = line.split()
            verts = list(map(int, parts[:-1]))  # derniers colonnes = tags parfois
            if v1 in verts and v2 in verts:
                n_removed += 1
                continue  # on supprime cet élément
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)

    with open(filename_out, "w") as f:
        f.writelines(out_lines)

    print(f"Suppression terminée : {n_removed} éléments supprimés.")
    print(f"Nouveau fichier : {filename_out}")


# Exemple d'utilisation
cut_link_in_mesh("plate.mesh", "plate_cut.mesh", v1=12, v2=34)
