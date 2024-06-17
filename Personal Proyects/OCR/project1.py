from PIL import Image
print("Librerias importadas con éxito")

im_file = r"C:\Users\alexe\Documentos\ProfessionalPortfolio\Personal Proyects\OCR\images\image1.jpg"

im = Image.open(im_file)
print(im) # Print image metadata
print(im.size) # Print image size

im.rotate(180).show() # Rotate the image but not resizing it
# im.show() # Show image

im.save(r"C:\Users\alexe\Documentos\ProfessionalPortfolio\Personal Proyects\OCR\images\results\img1.jpg")
print("Imagen guardada con éxito")