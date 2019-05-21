1. Create GlobusID
    * Go to [https://www.globusid.org/](https://www.globusid.org/).
    
2. Create ssh key
    * The workflow is described [here](https://support.opensciencegrid.org/support/solutions/articles/12000027675-generate-ssh-keys-and-activate-your-osg-login#step-1-generating-ssh-keys)
    * If you already have a file called id_rsa.pub, then don't overwrite it. Instead you can create another file.

3. Add ssh key to Globus
    * Login to [Globus](https://www.globusid.org/)
    * Click manage SSH and X.509 keys
    * Add the contents of the .pub file created above into the box.
    
4. Wait until OSG adds the key to the OSG connect website (will take hours).

5. 
