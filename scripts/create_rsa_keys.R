library(openssl)

# create and store a prive and public key for converting sample IDs -> anonymized/encryptied IDs
#
# sample IDs do not have personally identifiable information but just to be extra safe don't share
# info 
#
# Mark Stenglein Jan 26, 2022

# Generate keypair
key <- rsa_keygen()
pubkey <- as.list(key)$pubkey

# write out the private key
# this should not be shared and should have appropriate file permissions to remain unviewable by other users
write_pem(key, "rsa_key")

# write out the public key
# ok to share this / put on github
write_pem(pubkey, "rsa_key.pub")

