static void GetNewPrecipitateTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
                      int *nextAvailableTag)
{
        int     nextTag, thisDomain;

        nextTag    = *nextAvailableTag;
        thisDomain = home->myDomain;

/*
 *      If the old tag belonged to a different domain, we must
 *      give the precipitate a new tag.
 */
        if (oldTag->domainID != thisDomain) {

            for ( ; ; nextTag++) {

/*
 *              Extend the precipitatekeys array if necessary.
 */
                if (nextTag >= home->newPrecipitateKeyMax) {
                    ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_PRECIPITATEKEY_INC);
                }

                if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }
        } else {
/*
 *          The old tag belonged to this domain...  check if that tag
 *          is available for this run.
 *
 *          Extend the precipitatekeys array if necessary.
 */
            if (oldTag->index >= home->newPrecipitateKeyMax) {
                ExtendPrecipitateKeys(home, oldTag->index + NEW_PRECIPITATEKEY_INC);
            }

/*
 *          If the old tag is still available, use it and return
 *          to the caller.
 */
            if (home->precipitateKeys[oldTag->index] == (Precipitate_t *)NULL) {
                newTag->domainID = oldTag->domainID;
                newTag->index    = oldTag->index;
                if (newTag->index >= home->newPrecipitateKeyPtr) {
                    home->newPrecipitateKeyPtr = newTag->index + 1;
                }
                return;
            }

/*
 *          The old tag is no longer available, so just use
 *          the next available tag.
 */
            for ( ; ; nextTag++) {

/*
 *              Extend precipitate keys array if necessary
 */
                if (nextTag >= home->newPrecipitateKeyMax) {
                    ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_PRECIPITATEKEY_INC);
                }

                if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }
        }

        if (newTag->index >= home->newPrecipitateKeyPtr) {
            home->newPrecipitateKeyPtr = newTag->index + 1;
        }

        *nextAvailableTag = nextTag;

        return;
}
