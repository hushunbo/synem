import matplotlib.pyplot as plt
from skimage.filter import roberts, sobel


img_col = imread("brain.png")
edge_roberts = roberts(img)
edge_sobel = sobel(img)

#Normalize
edge_roberts = (edge_roberts - edge_roberts.min())/(edge_roberts.max() - edge_roberts.min())
edge_sobel = (edge_sobel - edge_sobel.min())/(edge_sobel.max() - edge_sobel.min())

#Clean
edge_roberts *= (edge_roberts>0.1) 
edge_sobel *= (edge_sobel>0.1)

# Red
brain_red = np.ones_like(img_col)
brain_red[:,:,1] = 1 - edge_sobel
brain_red[:,:,2] = 1 - edge_sobel
plt.figure()
imshow(brain_red)

# Black
brain_green= np.ones_like(img_col)
brain_green[:,:,0] = (1-edge_sobel)
brain_green[:,:,1] = (1-edge_sobel)
brain_green[:,:,2] = (1-edge_sobel)
plt.figure()
imshow(brain_green)

# Blue
brain_blue= np.ones_like(img_col)
brain_blue[:,:,0] = 1-edge_sobel
brain_blue[:,:,1] = 1-edge_sobel
plt.figure()
imshow(brain_blue)

# Negate
edge_roberts = edge_roberts.max() - edge_roberts
edge_sobel = edge_sobel.max() - edge_sobel




fig, (ax0, ax1) = plt.subplots(ncols=2)

ax0.imshow(edge_roberts, cmap=plt.cm.gray)
ax0.set_title('Roberts Edge Detection')
ax0.axis('off')

ax1.imshow(edge_sobel, cmap=plt.cm.gray)
ax1.set_title('Sobel Edge Detection')
ax1.axis('off')

plt.show()