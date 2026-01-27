import numpy as np
from typing import Union


def angle_between_vectors(
    v1: np.ndarray,
    v2: np.ndarray,
    degrees: bool = True
) -> float:
    """
    Calculate angle between two vectors.

    Args:
        v1: First 3D vector
        v2: Second 3D vector
        degrees: Return angle in degrees (default) or radians

    Returns:
        angle: Angle between vectors (0-180 degrees or 0-pi radians)
    """
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    cos_theta = np.dot(v1_norm, v2_norm)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    angle = np.arccos(cos_theta)

    if degrees:
        angle = np.degrees(angle)

    return float(angle)


def project_vector_onto_plane(
    vector: np.ndarray,
    plane_normal: np.ndarray
) -> np.ndarray:
    """
    Project vector onto plane.

    The projection is calculated as: projection = v - (v·n)n
    where n is the normalized plane normal.

    Args:
        vector: Vector to project
        plane_normal: Plane normal vector

    Returns:
        projected: Projected vector
    """
    normal_norm = plane_normal / np.linalg.norm(plane_normal)

    projection_length = np.dot(vector, normal_norm)

    projected = vector - projection_length * normal_norm

    return projected


def dihedral_angle(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    degrees: bool = True
) -> float:
    """
    Calculate dihedral angle defined by four points.

    The dihedral angle is the angle between two planes:
    - Plane 1: defined by points p0, p1, p2
    - Plane 2: defined by points p1, p2, p3

    Args:
        p0, p1, p2, p3: Four point coordinates
        degrees: Return angle in degrees (default) or radians

    Returns:
        angle: Dihedral angle (-180 to 180 degrees or -pi to pi radians)
    """
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    angle = np.arctan2(y, x)

    if degrees:
        angle = np.degrees(angle)

    return float(angle)


def signed_angle_between_vectors(
    v1: np.ndarray,
    v2: np.ndarray,
    normal: np.ndarray,
    degrees: bool = True
) -> float:
    """
    Calculate signed angle between two vectors in a plane.

    The sign is determined by the right-hand rule with respect to the normal vector.

    Args:
        v1: First vector
        v2: Second vector
        normal: Normal vector defining the plane orientation
        degrees: Return angle in degrees (default) or radians

    Returns:
        angle: Signed angle (-180 to 180 degrees or -pi to pi radians)
    """
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    normal_norm = normal / np.linalg.norm(normal)

    cos_theta = np.dot(v1_norm, v2_norm)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    cross = np.cross(v1_norm, v2_norm)
    sin_theta = np.dot(cross, normal_norm)

    angle = np.arctan2(sin_theta, cos_theta)

    if degrees:
        angle = np.degrees(angle)

    return float(angle)
