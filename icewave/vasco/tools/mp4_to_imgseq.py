import cv2
import os

#video_file = 'L:/Share_hublot/Data/0220/Drones/bernache/22-frac_001/DJI_20240220144844_0775_D.mp4'
#output_dir = 'C:/Users/Vasco Zanchi/Desktop/image_sequence_drone_20240220/'

def extract_imgseq_from_mp4(video_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    vidcap = cv2.VideoCapture(video_file)            # open the file
    if not vidcap.isOpened():
        print(f"cannot open {video_file}")
        return

    fps = vidcap.get(cv2.CAP_PROP_FPS)
    print("frame rate :", fps)

    count = 0
    while True:
        success, image = vidcap.read()               # read each frame
        if not success or image is None:            # end of stream
            break

        fname = os.path.join(output_dir, f"frame{count:04d}.jpg")
        cv2.imwrite(fname, image)

        if cv2.waitKey(1) == 27:                    # escape to quit
            break
        count += 1

    vidcap.release()
    print(f"extracted {count} frames into {output_dir}")


# === Paramètres ===
#input_folder = 'C:/Users/Vasco Zanchi/Desktop/image_sequence_drone_20240220_debut'     # à modifier
#output_folder = 'C:/Users/Vasco Zanchi/Desktop/images_decimees'    # à modifier
#num_frames_to_process = 700
#temporal_stride = 3  # Une image sur 3
#scale_factor = 0.1   # Réduction spatiale

# Créer le dossier de sortie s'il n'existe pas

def decimate_space_time_imgseq(input_folder,output_folder,num_frames_to_process,temporal_stride,scale_factor):

    os.makedirs(output_folder, exist_ok=True)

    # Traitement
    saved_frame_index = 0
    for i in range(0, num_frames_to_process, temporal_stride):
        filename = f"frame{i}.jpg"
        input_path = os.path.join(input_folder, filename)

        # Charger l'image
        img = cv2.imread(input_path)
        if img is None:
            print(f"Image non trouvée : {filename}")
            continue

        # Redimensionner l'image (diviser par 2 en largeur et hauteur)
        height, width = img.shape[:2]
        resized_img = cv2.resize(img, (width // 2, height // 2), interpolation=cv2.INTER_AREA)

        # Sauvegarder l'image décimée
        output_filename = f"frame{saved_frame_index}.jpg"
        output_path = os.path.join(output_folder, output_filename)
        cv2.imwrite(output_path, resized_img)

        saved_frame_index += 1

    print(f"Terminé. {saved_frame_index} images ont été sauvegardées dans {output_folder}.")
