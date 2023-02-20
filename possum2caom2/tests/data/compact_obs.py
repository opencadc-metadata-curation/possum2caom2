from glob import glob
from caom2pipe.manage_composable import read_obs_from_file, write_obs_to_file


def main():
    # make all the observations in this directory into one observation
    dir_listing = glob('./*.expected.xml')
    obs = None
    for entry in dir_listing:
        print(entry)
        temp_obs = read_obs_from_file(entry)
        if obs is None:
            obs = temp_obs
        else:
            for temp_plane in temp_obs.planes.values():
                if temp_plane.product_id in obs.planes.keys():
                    for temp_artifact in temp_plane.artifacts.values():
                        print(f'Adding artifact {temp_artifact.uri}')
                        obs.planes[temp_plane.product_id].artifacts.add(temp_artifact)
                else:
                    print(f'Adding plane {temp_plane.product_id}')
                    obs.planes.add(temp_plane)

    write_obs_to_file(obs, './compacted_obs.xml')


if __name__ == '__main__':
    main()
