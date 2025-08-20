import os
import numpy as np

def extract_value_with_indices(arr, layer=None, val_col='elev', mask_value=np.nan):
    import pandas as pd
    """
    Extract non-NaN values from array along with their indices.
    
    Args:
        arr: numpy array
        
    Returns:
        list: List of [index_tuple, value] pairs for non-NaN elements
    """
   
    # Get the indices of non-NaN elements
    indices = np.where(arr != mask_value)
    
    # Get the non-NaN values
    values = arr[arr != mask_value]
    
    # Combine indices and values into DataFrame format
    index_list = []
    elev_list = []
    
    for i in range(len(values)):
        if arr.ndim == 2:
            index_tuple = (indices[0][i], indices[1][i])
            if layer is not None:
                index_tuple = (layer, index_tuple[0], index_tuple[1])
        elif arr.ndim == 3:
            index_tuple = (indices[0][i], indices[1][i], indices[2][i])
        else:
            # Fallback for other dimensions, though not expected
            index_tuple = tuple(int(idx[i]) for idx in indices)
        
        index_list.append(tuple(map(int, index_tuple)))
        elev_list.append(float(values[i]))
    
    # Create DataFrame
    df = pd.DataFrame({
        'index': index_list,
        val_col: elev_list
    })

    return df


def get_interior_indices(arr, layer=None):
    import numpy as np
    from scipy.ndimage import binary_erosion
    """
    Extract boundary indices of an active domain (where values == 1).
    Boundary pixels are active pixels that have at least one inactive neighbor.
    
    Args:
        arr: numpy array where 1s represent active domain
        
    Returns:
        list: List of tuples representing boundary indices
    """
    # Create binary mask for active domain
    active_mask = (arr == 1)
    
    # Erode the mask - this removes boundary pixels
    eroded = binary_erosion(active_mask)
    
    # Boundary is active pixels minus eroded pixels
    boundary_mask = active_mask & ~eroded
    
    # Get indices where exterior mask is True
    result = get_indices(boundary_mask, layer)
    
    return boundary_mask, result

def get_exterior_indices(arr, layer=None):
    from scipy.ndimage import binary_dilation
    """
    Extract indices just outside the active domain boundary.
    These are inactive pixels that have at least one active neighbor.
    
    Args:
        arr: numpy array where 1s represent active domain
        
    Returns:
        list: List of tuples representing exterior indices
    """
    # Create binary mask for active domain
    active_mask = (arr == 1)
    
    # Dilate the mask - this expands the active area by 1 pixel
    dilated = binary_dilation(active_mask)
    
    # Exterior is dilated pixels minus original active pixels
    exterior_mask = dilated & ~active_mask.data
    
    # Get indices where exterior mask is True
    result = get_indices(exterior_mask, layer)
    
    return exterior_mask, result

def get_indices(arr, layer=None, value=False):
    # Get indices where exterior mask is True
    indices = np.where(arr)
    result = []
    for i in range(len(indices[0])):
        if arr.ndim == 2:
            index_tuple_2d = (indices[0][i], indices[1][i])
            ivalue = arr[index_tuple_2d]
            index_tuple = index_tuple_2d
            if layer is not None:
                index_tuple = (layer, index_tuple[0], index_tuple[1])
        elif arr.ndim == 3:
            index_tuple = (indices[0][i], indices[1][i], indices[2][i])
            ivalue = arr[index_tuple]
        else:
            # Fallback for other dimensions, though not expected
            index_tuple = tuple(int(idx[i]) for idx in indices)
            ivalue = arr[index_tuple]

        if value:
            result.append([tuple(map(int, index_tuple)), ivalue])
        else:
            result.append(tuple(map(int, index_tuple)))

    return result



def savedf2txt(df, filename):
    """
    Save a DataFrame to a text file with specified separator.
    
    Args:
        df: pandas DataFrame to save
        filename: Name of the output text file
        sim_ws: model workspace directory
    """
    df['k'] = df['index'].apply(lambda x: int(x[0] + 1))
    df['i'] = df['index'].apply(lambda x: int(x[1] + 1))
    df['j'] = df['index'].apply(lambda x: int(x[2] + 1))
    # delete df['index']  # remove index column if not needed
    df = df.drop(columns=['index'])
    df = df[['k', 'i', 'j'] + [col for col in df.columns if col not in ['k', 'i', 'j']]]

    df.to_csv(filename, sep=' ', header=False, index=False)

# plot layers
import matplotlib.pyplot as plt
import numpy as np

def plot_array_layers(arr, figsize=(15, 5), cmap='viridis', titles=None):
    """
    Plot each layer of a 3D array as subplots.
    
    Args:
        arr: 3D numpy array of shape (n_layers, height, width)
        figsize: Figure size (width, height)
        cmap: Colormap for the plots
        titles: List of titles for each subplot (optional)
    """
    n_layers = arr.shape[0]
    
    # Create subplots
    fig, axes = plt.subplots(1, n_layers, figsize=figsize)
    
    # Handle case where there's only one layer
    if n_layers == 1:
        axes = [axes]
    
    # Plot each layer
    for i in range(n_layers):
        im = axes[i].imshow(arr[i], cmap=cmap, vmin=np.nanmin(arr), vmax=np.nanmax(arr))
        
        # Set title
        if titles and i < len(titles):
            axes[i].set_title(titles[i])
        else:
            axes[i].set_title(f'Layer {i}')
        
    # Add colorbar
    plt.colorbar(im)
    
    plt.tight_layout()
    plt.show()